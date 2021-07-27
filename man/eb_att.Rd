\name{eb_att}
\alias{eb_att}

\title{
Estimate Energy Balancing Weights for the Average Treatment Effect on the Treated
}
\description{
\code{eb_att} is a wrapper for \code{\link{energybalance}} that focuses on using it to estimate energy balancing weights for the average treatment effect on the treated (ATT) with or more control groups. The energy distance between each non-focal treatment group and the focal treatment group is minimized.
}
\usage{
eb_att(sampleX, Z, sampleW = NULL, std = "studentized",
       focal = NULL, improved = FALSE, lambda = 0)
}
\arguments{
  \item{sampleX}{
a matrix or data frame of covariates composing the sample distribution to be weighted. See \code{\link{energybalance}} for how this is processed.
}
  \item{Z}{
a vector denoting treatment group membership.
}
  \item{sampleW}{
optional; a vector of sampling weights for the sample. The product of the estimated weights and \code{sampleW} will be the distribution-matching weights.
}
  \item{std}{
\code{character}; whether to standardize the covariates. If \code{"studentized"}, the distance matrix used will be the Euclidean distance with each variable scaled using its standard deviation in the focal treatment group distribution (weighted by \code{sampleW} if supplied). If \code{"mahalanobis"}, the distance matrix used will be the Mahalanobis distance computed using the covariance matrix of the focal treatment group distribution (weighted by \code{sampleW} if supplied). If \code{"none"}, the distance matrix used will be the Euclidean distance matrix. Default is \code{"studentized"}. Abbreviations allowed.
}
  \item{focal}{
a value of the \code{Z} to be considered the "focal" or "treated" group. If \code{NULL}, if \code{Z} is a factor, then \code{focal = levels(Z)[nlevels(Z)]}, and otherwise, \code{focal = max(Z)}.
}
  \item{improved}{
\code{logical}; whether to additionally balance the distributions between non-treatment groups when there are mor than two groups. If \code{TRUE}, the energy distance between non-focal treatment groups and between each non-focal treatment group and the focal treatment group distribution will be minimized. If \code{FALSE}, only the energy distance between each non-focal treatment group and the focal treatment group distribution will be minimized. Default is \code{FALSE} unlike other function in the package. With only two groups, this argument is ignored.
}
  \item{lambda}{
a penalty on the sum of the squared weights. \code{lambda/sum(Z!=focal)^2} times the sum of the squared weights is added to the objective function. Increasing \code{lambda} preserves the effective sample size of the weighted sample at the expense of balance.
}
}
\value{
A vector of weights, one for each row of \code{sampleX}. Weights for which \code{Z == focal} will be 1.
}

\author{
Noah Greifer
}

\seealso{
\code{\link{energybalance}} for a more flexible interface.

\code{\link{eb_ate}} for weights for estimating the average treatment effect (ATE).

\code{\link{eb_target}} for weights that match a sample to a target, which can also be used for estimating ATT weights.

\code{\link[WeightIt]{method_energy}} for the implementation of energy balancing in \emph{WeightIt}, which has additional options to restrict the difference in means between covariates.
}
\examples{

}
