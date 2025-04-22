import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse
from scipy import optimize
from scipy.sparse import lil_matrix, csr_matrix
from matplotlib import pyplot as plt
import pysam
from tqdm import tqdm

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

    ax.scatter(mu, var, c='k', label="observed RNA", rasterized=True, s=1)
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
