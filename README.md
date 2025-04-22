# scATAC-Seq Analysis Using a Poisson Model

- The motivation for this analysis was inspired by [this publication](https://www.nature.com/articles/s41592-023-02112-6) wherein it was observed that poisson models can be used with quantitative fragment counts to enhance the interpretability of scATAC-Seq data.

- While there is already a [tool](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/atac/PoissonVI.html) that exists within the scvi-tools suite to analyze such data, it is computationally expensive and inaccessible to researchers with lower-end computational power.

- Here, I am showcasing how orthogonal methods can be used to achieve similar results with standard single-cell tools that exist the scverse.

## Preprocessing

Publicly available data was first obtained from the 10X Genomics wesbsite:

- [Raw data](https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5)

- [Fragments file](https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz)

- [Fragments index](https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz.tbi)

We will also need an peak_annotation.tsv file as outlined [here](https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/analysis/peak-annotations).

In the misc_scripts folder, there is home-brew command line tool that can generate this file by taking in two fiiles:

1. gtf file containing transcriptomic coordinates
2. filtered .h5 matrix output from cellRanger

Below is an example of its usage from a jupyter notebook:

<pre> ```python ~/misc_scripts/map_peaks_to_genes.py \
  --input_file atac_v2_pbmc_10k_filtered_peak_bc_matrix.h5 \ 
  --gtf_file ~/refdata-gex-GRCh38-2020-A/genes/genes.gtf </pre>

