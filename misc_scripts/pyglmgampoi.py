import numpy as np
import rpy2
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
from rpy2.robjects import Formula
from rpy2.robjects import IntVector
from rpy2.robjects import pandas2ri
from rpy2.robjects import r
from rpy2.robjects.packages import importr

pandas2ri.activate()

def fit_glmgp_offset(y, coldata, log_umi, overdispersion, design="~ 1"):
    glmgp = importr("glmGamPoi")
    y_ro = np.asmatrix(y)
    fit = glmgp.glm_gp(
        data=y_ro,
        design=Formula(design),
        col_data=coldata,
        size_factors=False,
        offset=log_umi,
        overdispersion=overdispersion
    )
    overdispersions = fit.rx2("overdispersions")
    mu = fit.rx2("Mu")
    beta = fit.rx2("Beta")
    return {
        "overdispersions": np.squeeze(overdispersions),
        "mu": np.squeeze(mu.mean(1)),
        "theta" : np.vstack((1 / overdispersions, mu.mean(1) / 1e-4)).min(axis=0),
        "Intercept" : np.squeeze(beta),
        "log10_umi" : np.log(10)
    }

def bw_SJr(y, bw_adjust=3):
    if rpy2 is None:
        raise ImportError("bw_SJr requires rpy2 which is not installed.")
    rpy2.robjects.numpy2ri.activate()
    stats = importr("stats")
    base = importr("base")
    return np.asarray(stats.bw_SJ(y)) * bw_adjust

def ksmooth(genes_log_gmean, genes_log_gmean_step1, col_to_smooth, bw):
    if rpy2 is None:
        raise ImportError("ksmooth requires rpy2 which is not installed.")
    rpy2.robjects.numpy2ri.activate()
    stats = importr("stats")
    base = importr("base")
    x_points = base.pmax(genes_log_gmean, base.min(genes_log_gmean_step1))
    x_points = base.pmin(x_points, base.max(genes_log_gmean_step1))
    o = base.order(x_points)
    dispersion_par = stats.ksmooth(
        x=genes_log_gmean_step1,
        y=col_to_smooth,
        x_points=x_points,
        bandwidth=bw,
        kernel="normal",
    )
    dispersion_par = dispersion_par[dispersion_par.names.index("y")]
    return {"smoothed": dispersion_par, "order": o}