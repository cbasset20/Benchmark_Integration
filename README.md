# Benchmark_Integration
Compare single cell integration strategies. Based on the [scib](https://github.com/theislab/scib) package. See [this tutorial](https://scib-metrics.readthedocs.io/en/stable/notebooks/lung_example.html). The main challenge is in converting seurat objects into anndata objects so they can be compared using scib_metrics. To solve this, the seurat object must first be deconstructed into its component parts. These are then read in via python to reconstruct an anndata object. See [this tutorial](https://jr-leary7.github.io/quarto-site/tutorials/Seurat_AnnData_Conversion.html) for additional information.

This guide provides a systematic workflow for benchmarking scRNA-seq integration algorithms. It resolves R-to-Python interoperability by extracting the underlying count matrices, metadata, and reduction embeddings from Seurat. This ensures the accurate reconstruction of UMAPs and integrated spaces within a native AnnData object, seamlessly bridging the gap to Python-based evaluation pipelines.

#Installation
# in here I will write step by step how i installed my dependencies and setup my environments.
#before running python installation

#https://visualstudio.microsoft.com/visual-cpp-build-tools/

#Download and run the installer

#In the installer: ✅ Select “C++ build tools” ✅ Make sure these are checked:

#MSVC v14.x (x64/x86)
#Windows 10 or 11 SDK
#C++ CMake tools (optional but recommended)
#restart computer after installation


Python
```bash
conda create --name scib-pipeline-R4.0 python=3.9
conda activate scib-pipeline-R4.0
pip install jupyter==1.1.1
pip install faiss-cpu 
pip install swig==4.4.1  # see: https://github.com/facebookresearch/faiss/issues/4481
pip install leidenalg==0.10.2  # scib dependency issue
pip install scib==1.1.7  # note: this fails to build under python=3.10 due to llvmlite
pip install scib-metrics==0.5.1
pip install scanorama==1.7.4
pip install pyliger==0.2.4
pip install harmony-pytorch==0.1.8
pip install scvi-tools==1.1.6
pip install anndata==0.10.8  # 0.10.9 causes import errors
```

R
```bash
# Azimuth install for Seurat v5
BiocManager::install("CNEr")  # TFBSTools dependency
devtools::install_github("ge11232002/TFMPvalue")  # TFBSTools dependency #NOTE: NEEDS R version >= 4.5.0 #install R 4.5.2 and Rtools45
BiocManager::install("Biostrings", force=TRUE)  # TFBSTools dependency
BiocManager::install("TFBSTools", force = TRUE)
BiocManager::install("DirichletMultinomial", force=TRUE)  # sudo apt install libgsl-dev
install.packages("hdf5r", configure.args = "--with-hdf5=/usr/bin/h5cc")  # see: https://github.com/hhoeflin/hdf5r/issues/112#issuecomment-453349767
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")  # seurat-disk dependency
devtools::install_github('satijalab/seurat-data')
remotes::install_github("mojaveazure/seurat-disk")
BiocManager::install(version = '3.22')
BiocManager::install("EnsDb.Hsapiens.v86", version = "3.22")#azimuth dependency
remotes::install_github('satijalab/azimuth', ref = 'master')
remotes::install_github("bnprks/BPCells/r")  # see below for mac
BiocManager::install("biomaRt", force=TRUE, lib="")
install.packages("this.path")  # don't use for Windows
install.packages("qs2")
install.packages("optparse")
install.packages("logr")
install.packages("import")
```

#To follow this tutorial:https://jr-leary7.github.io/quarto-site/tutorials/Seurat_AnnData_Conversion.html
#install these packages:
 BiocManager::install("scRNAseq")
install.packages("dplyr")
 install.packages("paletteer")
remotes::install_github("jessegmeyer/scLANE")
