\name{ShowResiduals}
\alias{ShowResiduals}
\title{Plots MOC Z-score residuals.}
\usage{
ShowResiduals(CoExpr)
}
\description{
Visualize the residuals of the MOC Z-score vs. the MOC Ratio.}
\arguments{
\item{CoExpr}{\code{CoExpress()} output \code{data.frame}.}
}
\value{
Returns a grob showing the linear regression of MOC Z ~ MOC Ratio, and if calculated, the adjusted Z-scores.
}
\examples{
# ShowResiduals(CoExpr)
}
