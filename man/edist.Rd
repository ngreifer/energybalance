\name{edist}
\alias{edist}
\title{
Weighted Energy Distance
}
\description{
Computes the (weighted) energy distance between two distributions.
}
\usage{
edist(sampleX, targetX = NULL, sampleW = NULL, targetW = NULL,
      Z = NULL, std = "studentized")
}
\arguments{
  \item{sampleX}{
a matrix of covariates composing the sample distribution.
}
  \item{targetX}{
a matrix of covariates composing the target distribution. If \code{NULL}, \code{sampleX} will be used.
}
  \item{sampleW}{
optional; a vector of weights for the sample.
}
  \item{targetW}{
optional; a vector of weights for the target.
}
  \item{Z}{
a vector denoting group membership. See Details.
}
  \item{std}{
\code{character}; whether to standardize the covariates. If \code{"studentized"}, the distance matrix used will be the Euclidean distance with each variable scaled using its standard deviation in the target distribution (weighted by \code{targetW} if supplied). If \code{"mahalanobis"}, the distance matrix used will be the Mahalanobis distance computed using the covariance matrix of the target distribution (weighted by \code{targetW} if supplied). If \code{"none"}, the distance matrix used will be the Euclidean distance matrix. Default is \code{"studentized"}. Abbreviations allowed.
}
}
\details{
\code{edist()} computes the (weighted) energy distance between a sample and a target distribution or between two samples distributions. When \code{Z} is \code{NULL}, the (weighted) energy distance is computed between the sample represented by \code{sampleX} and the target represented by \code{targetX}. If \code{std} is \code{TRUE}, the standardization occurs using \code{targetX} (and \code{targetW} if supplied). When \code{Z} is not \code{NULL}, the (weighted) energy distance is computed between the groups denoted by \code{Z}. If \code{std} is \code{TRUE}, the standardization occurs using \code{targetX} and \code{targetW}, but otherwise, they are ignored.

To compute the energy distance between a weighted sample and a target distribution, use
\preformatted{edist(sampleX, targetX = targetX, sampleW = sampleW)}

To compute the energy distance between two weighted groups denoted by \code{Z}, use
\preformatted{edist(sampleX, sampleW = sampleW, Z = Z)}

To compute the criterion minimized by \code{eb_ate()} with \code{improved = TRUE}, use
\preformatted{
edist(sampleX[Z==1,], targetX = sampleX, sampleW = sampleW[Z==1]) +
  edist(sampleX[Z==0,], targetX = sampleX, sampleW = sampleW[Z==0]) +
  edist(sampleX, sampleW = sampleW, Z = Z)}
}
\value{
A single value denoting the (weighted) energy distance between the two distributions.
}
\author{
Noah Greifer
}
\seealso{
\code{\link{energybalance}} to estimate weights that minimize the energy distance.
}
\examples{

}
