import scanpy as sc
import numpy as np
import pandas as pd
import numba as nb
import seaborn as sns
import scipy.sparse
from scipy import optimize
from scipy.sparse import lil_matrix, csr_matrix
from matplotlib import pyplot as plt
import pysam
from tqdm import tqdm
from scipy.stats import norm
from scipy.special import expit
from statsmodels.stats.multitest import multipletests
from sklearn.exceptions import ConvergenceWarning
from sklearn import metrics
from sklearn.metrics import confusion_matrix, classification_report
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import train_test_split, permutation_test_score
from sklearn.preprocessing import StandardScaler
from sklearn.inspection import permutation_importance
from numpy.typing import NDArray

def residual_vars(adata, clip=None, theta=np.inf, layer=None, chunk_size=10000):
    
    #get depth-adjusted expected counts without zero-inflation
    counts = sc.get._get_obs_rep(adata, layer=layer)

    sum_genes = np.array(counts.sum(axis=0)).squeeze()
    sum_cells = np.array(counts.sum(axis=1)).squeeze()
    sum_total = np.sum(sum_genes).squeeze()
    mu = np.outer(sum_cells, sum_genes) / sum_total

    # Allocate storage for variance
    variances = np.zeros(counts.shape[1], dtype=np.float32)
    
    for start in range(0, counts.shape[1], chunk_size):
        end = min(start + chunk_size, counts.shape[1])
        print(f"Processing peaks {start} to {end - 1}")

        counts_chunk = counts[:, start:end]
        mu_chunk = mu[:, start:end]

        variance_chunk = mu_chunk + np.divide(mu_chunk**2, theta)
        res_chunk = np.divide(counts_chunk - mu_chunk, np.sqrt(variance_chunk))

    # prepare clipping
        clip_val = np.sqrt(mu.shape[0]) if clip is None else clip
        res_chunk = np.clip(res_chunk, a_min=-clip_val, a_max=clip_val)

        # Compute variance in chunks
        variances[start:end] = np.var(res_chunk, axis=0)
        
    return pd.DataFrame(variances, index=adata.var_names.tolist(), columns=['residual_variance'])



def make_fragment_matrix (adata, peaks, fragments_file):
    
    fragments = pysam.TabixFile(fragments_file, parser=pysam.asBed()) #use pysam to create a tabix file
    
    cell_idx = {cell: i for i, cell in enumerate(adata.obs.index.unique())} #get the cell index
    f_idx = {f: i for i, f in enumerate(adata.var.index)} #get the fragment index
    
    rows, cols, values = [], [], [] #empty lists to store the rows, columns and values of the sparse matrix

    n_cells = len(cell_idx)
    n_peaks = len(f_idx)

    fragment_matrix = lil_matrix((n_cells, n_peaks), dtype=np.float32) #initialize the sparse matrix
    
    for i in tqdm(range(peaks.shape[0])):
        f = peaks.iloc[i] #process each individual peak
        fr = fragments.fetch(f.Chromosome, f.Start, f.End) #fetch fragments in the peak using pysam tabix file
        df = pd.DataFrame(
                [(x.contig, x.start, x.end, x.name) for x in fr],
                columns=["Chromosome", "Start", "End", "Cell"],
            )

        feature = f.Chromosome + ":" + str(f.Start) + "-" + str(f.End) #needs to be in the format of Chromosome:Start-End
        col_index = f_idx.get(feature) #fetch the index of the feature

        c_dict = df.Cell.value_counts().to_dict() #count the number of fragments in each cell

        for cell, count in c_dict.items():
            row_index = cell_idx.get(cell) #fetch the index of the cell
            if row_index is not None:
                fragment_matrix[row_index, col_index] = count #add the count to the sparse matrix

    # Convert to CSR format for efficient computation
    fragment_matrix = fragment_matrix.tocsr()

    return fragment_matrix

def mean_var_curvefit (adata, mu, var, nb=False, ax=None, title=None):

    if ax is None:
        fig, ax = plt.subplots()

    phi_hat, _ = optimize.curve_fit(lambda m, phi: m + phi * m ** 2, mu.values , var.values) #fitting mean variance relationship to negative binomial

    mm = np.logspace(np.log10(mu.values.min()), np.log10(mu.values.max()), num=len(mu))

    ax.scatter(mu, var, c='k', label="observed fragments", rasterized=True, s=1)
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.plot(mm, mm, label = 'Poisson', c='r') #poisson
    if nb:
        ax.plot(mm, (mm + phi_hat * mm ** 2), label = "Negative Binomial") #negative binomial


    ax.legend()
    if title is None:
        ax.set_title('Mean-Variance relationship')
    else:
        ax.set_title(title)

     # If we created the axis, show the plot
    if ax is None:
        plt.tight_layout()
        plt.show()

def depth_adjusted_mean_var(adata, layer=None, chunk_size=1000):

    #Step 1: #get depth-adjusted expected counts without zero-inflation

    counts = sc.get._get_obs_rep(adata, layer=layer) #extract counts matrix

    sum_genes = np.array(counts.sum(axis=0)).squeeze() #sum of all genes
    sum_cells = np.array(counts.sum(axis=1)).squeeze() #sum of all cells

    sum_total = np.sum(sum_genes).squeeze() #sum of all counts
    mu = np.outer(sum_cells, sum_genes) / sum_total #use linear algebra to calculate depth-adjusted means

    variances = np.zeros(counts.shape[1], dtype=np.float32) #allocate storage for variances

    #Step 2: Calculate variances in chunks to save memory

    for start in range(0, counts.shape[1], chunk_size):
        end = min(start + chunk_size, counts.shape[1])
        print(f"Processing genes {start} to {end - 1}")

        counts_chunk = counts[:, start:end]
        mu_chunk = mu[:, start:end]

        variances[start:end] = np.var(counts_chunk - mu_chunk, axis=0)

    #Step 3: Return results as a dictionary

    return {
        'gene_names' : adata.var_names,
        'variances': variances.reshape(1,-1),
        "da_means": mu
    }


def NB_drop_prob(mu, phi):
    #calculate the drop-out probability under a negative binomial model
    #if phi = 0, this reduces to a Poisson model
        if phi == 0:
            return np.exp(-mu)
        return ((1. / phi) / (mu + (1. / phi))) ** (1. / phi)

def lr_classifier(adata, hvgs_list, layer, cluster_key, solver, max_iter, train_size, plot=True):
    
    subset = adata[:, adata.var_names.isin(hvgs_list)].copy()
    counts = sc.get._get_obs_rep(subset, layer=layer)
    genes = subset.var_names
    key = subset.obs[cluster_key]
    print(f'Finished subsetting by genes')
    
    print('splitting into training and test datasets')
    X_train, X_test, y_train, y_test = train_test_split(counts, key, train_size=train_size, random_state=42)
    scaler = StandardScaler() # zero mean and unit variance 
    X_train = scaler.fit_transform(X_train.toarray()) #scale training data
    X_train = np.clip(X_train, -np.sqrt(X_train.shape[0]), np.sqrt(X_train.shape[0])) #clip extreme values
    print(f'Finished normalizing training data')

    #init classifier
    classifier = OneVsRestClassifier(LogisticRegression(solver = solver,
                                    C = 1,
                                    n_jobs = -1,
                                    random_state = 42))
    
    classifier.fit(X_train, y_train) #fit to training data
    classifier.features = genes #set features as genes
    print(f'Finished logistic regression fitting on training data')
    
    means_ = scaler.mean_ #extract mean
    sd_ = scaler.scale_ #extract sd
    X_test = (X_test.toarray() - means_) / sd_ #scale test data
    X_test = np.clip(X_test, -np.sqrt(X_test.shape[0]), np.sqrt(X_test.shape[0])) #clipping again
        
    scores = classifier.decision_function(X_test) #extract scores
    y_prob = expit(scores) #convert to probabilities
    
    if plot:
        fig, axes = plt.subplots(1, 2, figsize=(6, 3))
        for i, cell_type in enumerate(classifier.classes_):
            fpr, tpr, _ = metrics.roc_curve(y_test == cell_type, y_prob[:, i])
            axes[0].plot(fpr, tpr, lw=2)

        axes[0].plot([0, 1], [0, 1], color='k', ls=':')
        axes[0].set_xlabel('FPR')
        axes[0].set_ylabel('TPR')
        
        y_pred = classifier.classes_[scores.argmax(axis=1)]
        cmat = confusion_matrix(y_test, y_pred)
        cmat = cmat/cmat.sum(1)
        sns.heatmap(cmat, cmap='viridis', ax=axes[1])
        print(classification_report(y_test, y_pred))
    
    # return scaler, classifier, counts, key

def lr_classifier(adata, hvgs_list, layer, cluster_key, solver, max_iter, train_size, plot=True):
    
    subset = adata[:, adata.var_names.isin(hvgs_list)].copy()
    counts = sc.get._get_obs_rep(subset, layer=layer)
    genes = subset.var_names
    key = subset.obs[cluster_key]
    print(f'Finished subsetting by genes')
    
    print('splitting into training and test datasets')
    X_train, X_test, y_train, y_test = train_test_split(counts, key, train_size=train_size, random_state=42)
    scaler = StandardScaler() # zero mean and unit variance 
    X_train = scaler.fit_transform(X_train.toarray()) #scale training data
    X_train = np.clip(X_train, -np.sqrt(X_train.shape[0]), np.sqrt(X_train.shape[0])) #clip extreme values
    print(f'Finished normalizing training data')

    #init classifier
    classifier = OneVsRestClassifier(LogisticRegression(solver = solver,
                                    C = 1,
                                    n_jobs = -1,
                                    random_state = 42))
    
    classifier.fit(X_train, y_train) #fit to training data
    classifier.features = genes #set features as genes
    print(f'Finished logistic regression fitting on training data')
    
    means_ = scaler.mean_ #extract mean
    sd_ = scaler.scale_ #extract sd
    X_test = (X_test.toarray() - means_) / sd_ #scale test data
    X_test = np.clip(X_test, -np.sqrt(X_test.shape[0]), np.sqrt(X_test.shape[0])) #clipping again
        
    scores = classifier.decision_function(X_test) #extract scores
    y_prob = expit(scores) #convert to probabilities
    
    if plot:
        fig, axes = plt.subplots(1, 2, figsize=(6, 3))
        for i, cell_type in enumerate(classifier.classes_):
            fpr, tpr, _ = metrics.roc_curve(y_test == cell_type, y_prob[:, i])
            axes[0].plot(fpr, tpr, lw=2)

        axes[0].plot([0, 1], [0, 1], color='k', ls=':')
        axes[0].set_xlabel('FPR')
        axes[0].set_ylabel('TPR')
        
        y_pred = classifier.classes_[scores.argmax(axis=1)]
        cmat = confusion_matrix(y_test, y_pred)
        cmat = cmat/cmat.sum(1)
        sns.heatmap(cmat, cmap='viridis', ax=axes[1])
        print(classification_report(y_test, y_pred))

