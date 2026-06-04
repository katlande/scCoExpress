\name{plotCoExpr}
\alias{plotCoExpr}
\title{Plots Co-expression clusters.}
\usage{
plotCoExpr(co, label_df=NULL, fontsize=7, lineweight=0.5)
}
\description{Visualize co-expression clusters.}
\arguments{
\item{co}{CoExpress function output data.frame.}
\item{label_df}{Optional data.frame with three columns: gene=1, id=2, colour=3. Allows users to colour code axis text.}
\item{fontsize}{Numeric, axis text size; auto=7.}
\item{lineweight}{Numeric, line thickness; auto=0.5.}
}
\value{
Returns a grob, showing a clustered co-expression matrix of all the gene pairs tabulated in by \code{CoExpress()}.
}
\examples{
# plotCoExpr(CoExpr)
}
