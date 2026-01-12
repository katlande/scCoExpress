\name{plotCoExpr}
\alias{plotCoExpr}
\title{Plots Co-expression clusters.}
\usage{
plotCoExpr(co)
}
\description{
Visualize co-expression clusters.}
\arguments{
\item{co}{\code{CoExpress()} output \code{data.frame}.}
\item{label_df}{\code{Optional \code{data.frame} with three columns: gene (1), id (2), colour (3). Allows users to colour code axis text.}
\item{fontsize}{\code{numeric}, axis text size; auto=7.}
\item{linewight}{\code{numeric}, line thickness; auto=0.5.}
}
\value{
Returns a grob, showing a clustered co-expression matrix of all the gene pairs tabulated in by \code{CoExpress()}.
}
\examples{
# plotCoExpr(CoExpr)
}
