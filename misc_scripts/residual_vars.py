import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse
from scipy.sparse import lil_matrix, csr_matrix

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