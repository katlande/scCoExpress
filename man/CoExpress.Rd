\name{CoExpress}
\alias{CoExpress}
\title{Get Co-expression information from a list of genes.}
\usage{
CoExpress(obj, target_genes, gene2=NULL, seuratAssay=NULL, seuratSlot="data", nPartitions=10, nPermutations=50, BkgdGeneExpr=NULL, topExcl=0.98, bottomExcl=0.005, local.perms=6, skip.extremes=T, seurat=5)
}
\description{
Get cell-level co-expression values on a set of gene pairs from a single cell, single nuclei, or spatial dataset..}
\arguments{

\item{obj}{Seurat Object}
\item{target_genes}{Vector of gene names to query.}
\item{gene2}{Target gene if comparing \code{target_genes} against a single gene only, or NULL if querying all pairwise comparisons of \code{target_genes}. Default=NULL}
\item{seuratAssay}{Seurat Assay to use; if NULL, uses defaultAssay.}
\item{seuratSlot}{Seurat assay slot to use, defaults to "data."}
\item{nPartitions}{Number of gene partitions to use for background querying. If set to NULL, uses local backgrounds for all genes.}
\item{nPermutations}{Number of permutations for parition backgrounds, or for local backgrounds if nParitions is NULL.}
\item{BkgdGeneExpr}{getGeneExpr() output object, runs internally if NULL.}
\item{topExcl}{Value between 0-1, percentile above which cells will be set as "extreme" genes.}
\item{bottomExcl}{Value between 0-1, minimum fraction of cells a gene needs to be expressed in, below which genes will be set as "extreme."}
\item{local.perms}{In partition mode if \code{skip.extremes==FALSE}, the number of permutations to run on high and low expressed genes locally}
\item{skip.extremes}{\code{Boolean}, whether or not to skip tabulating co-expression on comparisons with very highly or very lowly expressed genes. If \code{TRUE}, runs extreme genes in local mode.}
\item{seurat}{Top level Seurat version to use as an \code{integer}. E.g., for Seurat 5.1.0, v=5.}
}
\value{
Returns Co-Expression results from a specified set of genes in a Seurat Object.
}
\examples{
# Local backgrounds for Apoe, Cd44, and Ltb vs. Ptprc, comparisons of extreme genes are not tabulated:
# CoEx_1 <- CoExpress(SeuratObj, c("Cd44", "Ltb", "Apoe"), gene2 = "Ptprc", nPermutations=10, nPartitions = NULL)

# Backgrounds for Apoe, Cd44, and Ltb vs. Ptprc using 7 partitions, comparisons of extreme genes are not tabulated:
# CoEx_2 <- CoExpress(SeuratObj, c("Cd44", "Ltb", "Apoe"), gene2 = "Ptprc", nPermutations=50, nPartitions = 7)

# Backgrounds for all pariwise comparisons of Apoe, Cd44, and Ltb using 10 partitions, comparisons of extreme genes use local backgrounds:
# CoEx_3 <- CoExpress(SeuratObj, c("Cd44", "Ltb", "Apoe"), nPermutations=50, nPartitions = 10, skip.extremes = F)
}
