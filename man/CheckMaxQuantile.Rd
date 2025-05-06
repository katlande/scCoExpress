\name{CheckMaxQuantile}
\alias{CheckMaxQuantile}
\title{Examine abundance of genes in your seuratObj at different.}
\usage{
CheckMaxQuantile(exprDF, quantiles=c(0.9, 0.96, 0.98))
}
\description{
Visualize the distribution of \%cells expressing the genes in a \code{SeuratObject}.}
\arguments{
\item{exprDF}{\code{getGeneExpr()} output \code{data.frame}.}
\item{quantiles}{Vector of values from 0-1, only uses the first 3. Percentiles to highlight on output plot.}
}
\value{
Returns a grob showing the density of genes by \% cells expressed.
}
\examples{
# CheckMaxQuantile(exprDF)
}
