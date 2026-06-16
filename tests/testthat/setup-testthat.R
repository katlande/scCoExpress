devtools::install_github('satijalab/seurat-data')
library(SeuratData)
library(Seurat)

# seurat object to use for testing:
InstallData('pbmc3k')
pbmc <- LoadData('pbmc3k')

# Some semi-random genes to use for the tests:
GOIs <- c("TSR2", "MPC2", "FCN1", "GPX4", "CD74", "HLA-DQB1")
