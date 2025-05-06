\name{getGeneExpr}
\alias{getGeneExpr}
\title{Gets gene abundances.}
\usage{
getGeneExpr(SeuratObj, assay=NULL, slot="data", v=5)
}
\description{
Get relative gene expression abundances from a \code{SeuratObject}.}
\arguments{
\item{SeuratObj}{A \code{SeuratObject}}
\item{assay}{Assay of the \code{SeuratObject} to use. If NULL, uses the DefaultAssay of the Seurat object.}
\item{slot}{Slot of \code{assay} to use; assumes "data" unless otherwise specified.}
\item{v}{Top level Seurat version to use as an \code{integer}. E.g., for Seurat 5.1.0, v=5.}
}
\value{
Returns a \code{data.frame} with the fraction of cells expressing each gene in a seurat object
}
\examples{
# gene_info <- getGeneExpr(SeuratObj)
}
