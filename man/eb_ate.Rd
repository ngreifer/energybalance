\name{eb_ate}
\alias{eb_ate}

\title{
Estimate Energy Balancing Weights for the Average Treatment Effect
}
\description{
\code{eb_ate} is a wrapper for \code{\link{energybalance}} that focuses on using it to estimate energy balancing weights for the average treatment effect (ATE). The energy distance between each treatment and the full sample is minimized.
}
\usage{
eb_ate(sampleX, Z, sampleW = NULL, std = "studentized",
       improved = TRUE)
}
\arguments{
  \item{sampleX}{
a matrix of covariates composing the sample distribution to be weighted.
}
  \item{Z}{
a vector denoting treatment group membership.
}
  \item{sampleW}{
optional; a vector of sampling weights for the sample. The product of the estimated weights and \code{sampleW} will be the distribution-matching weights.
}
  \item{std}{
\code{character}; whether to standardize the covariates. If \code{"studentized"}, the distance matrix used will be the Euclidean distance with each variable scaled using its standard deviation in the target distribution (weighted by \code{sampleW} if supplied). If \code{"mahalanobis"}, the distance matrix used will be the Mahalanobis distance computed using the covariance matrix of the target distribution (weighted by \code{sampleW} if supplied). If \code{"none"}, the distance matrix used will be the Euclidean distance matrix. Default is \code{"studentized"}. Abbreviations allowed.
}
  \item{improved}{
\code{logical}; whether to additional balance the distributions between treatment groups. If \code{TRUE}, the energy distance between treatment groups and between each treatment group and the full sample distribution will be minimized. If \code{FALSE}, only the energy distance between each treatment group and the full sample distribution will be minimized. Default is \code{TRUE} as recommended by Huling and Mak (2020).
}
}
\value{
A vector of weights, one for each row of \code{sampleX}.
}

\author{
Noah Greifer
}

\seealso{
\code{\link{energybalance}} for a more flexible interface that allows matching both groups to a target distribution other than the full sample.

\code{\link[WeightIt]{method_energy}} for the implementation of energy balancing in \emph{WeightIt}, which has additional options to restrict the difference in means between covariates and can incorproate a penalty on the dispersion of the weights.
}
\examples{

}