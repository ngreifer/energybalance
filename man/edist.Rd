\name{edist}
\alias{edist}
\alias{edist_ate}
\title{
Weighted Energy Distance
}
\description{
\code{edist()} computes the (weighted) energy distance between two distributions. \code{edist_ate()} is a wrapper for \code{edist()} specifically for computing the energy distance minimized when using \code{\link{eb_ate}} to estimate weights for the ATE.
}
\usage{
edist(sampleX, targetX = NULL, sampleW = NULL, targetW = NULL,
      Z = NULL, std = "studentized")

edist_ate(sampleX, Z, sampleW = NULL, targetW = NULL,
          std = "studentized", improved = TRUE)
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
a vector denoting treatment group membership. See Details.
}
  \item{std}{
\code{character}; whether to standardize the covariates. If \code{"studentized"}, the distance matrix used will be the Euclidean distance with each variable scaled using its standard deviation in the target distribution (weighted by \code{targetW} if supplied). If \code{"mahalanobis"}, the distance matrix used will be the Mahalanobis distance computed using the covariance matrix of the target distribution (weighted by \code{targetW} if supplied). If \code{"none"}, the distance matrix used will be the Euclidean distance matrix. Default is \code{"studentized"}. Abbreviations allowed.
}
  \item{improved}{
\code{logical}; whether to include the pairwise energy distances between treatment groups. If \code{TRUE}, the combined energy distance between treatment groups and between each treatment group and the full sample distribution will be computed. If \code{FALSE}, only the combined energy distance between each treatment group and the full sample distribution will be computed. Default is \code{TRUE}.
}
}
\details{
\code{edist()} computes the (weighted) energy distance between a sample and a target distribution or between two samples distributions. When \code{Z} is \code{NULL}, the (weighted) energy distance is computed between the sample represented by \code{sampleX} and the target represented by \code{targetX}. If \code{std} is \code{TRUE}, the standardization occurs using \code{targetX} (and \code{targetW} if supplied). When \code{Z} is not \code{NULL}, the (weighted) energy distance is computed between the groups denoted by \code{Z}. If \code{std} is \code{TRUE}, the standardization occurs using \code{targetX} and \code{targetW}, but otherwise, they are ignored.

To compute the energy distance between a weighted sample and a target distribution, use
\preformatted{edist(sampleX, targetX = targetX, sampleW = sampleW)}

To compute the energy distance between two weighted groups denoted by \code{Z}, use
\preformatted{edist(sampleX, sampleW = sampleW, Z = Z)}

\code{edist_ate()} computes the combined (weighted) energy distance between each treatment group and the full sample and, when \code{improved = TRUE}, between each pair of treatment groups. This is the quantity minimized when using \code{\link{eb_ate}} and is the sum of several energy distance components.
\preformatted{
edist_ate(sampleX, Z = Z, sampleW = sampleW,
          improved = TRUE)
}
is equivalent to
\preformatted{
edist(sampleX, sampleW = sampleW, Z = Z) +
  edist(sampleX[Z==1,], targetX = sampleX,
        sampleW = sampleW[Z==1]) +
  edist(sampleX[Z==0,], targetX = sampleX,
        sampleW = sampleW[Z==0])
  }
when the tretament variable is binary.
}
\value{
A single value denoting the (weighted) energy distance between the two distributions. For \code{edist_ate()}, the output will have the attribute \code{"components"} containing all the component energy distances summed to compute the returned value.
}
\author{
Noah Greifer
}
\seealso{
\code{\link{energybalance}} to estimate weights that minimize the energy distance, \code{\link{eb_ate}} to estimate weights that minimize the combined energy distance for estimating the ATE.
}
\examples{

}
