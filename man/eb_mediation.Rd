\name{eb_mediation}
\alias{eb_mediation}
\title{
Estimate Energy Balancing Weights for Mediation Analysis
}
\description{
\code{eb_mediation} is a wrapper for \code{\link{energybalance}} that focuses on using it to estimate energy balancing weights for the mean cross-world potential outcome in mediation analysis. Estimation proceeds in two steps. First, the control group covariate distribution is weighted to resemble the full sample. Second, the treated group covariate and mediator distribution is weighted to resemble the weighted control group. Weights are produced for the treated group.
}
\usage{
eb_mediation(sampleX, Z, M, sampleW = NULL,
             std = "studentized", lambda = 0)
}
\arguments{
  \item{sampleX}{
a matrix of covariates composing the sample distribution to be weighted.
}
  \item{Z}{
a vector denoting treatment group membership. The group with \code{Z==1} will be the treated group for which the weights are being estimated, and the group with \code{Z==0} will be the control group that the weighted treatment group is meant to resemble.
}
  \item{M}{
a vector containing the values of the mediator.
}
  \item{sampleW}{
optional; a vector of sampling weights for the sample. The product of the estimated weights and \code{sampleW} will be the distribution-matching weights.
}
  \item{std}{
\code{character}; whether to standardize the covariates. If \code{"studentized"}, the distance matrix used will be the Euclidean distance with each variable scaled using its standard deviation in the target distribution (weighted by \code{targetW} if supplied). If \code{"mahalanobis"}, the distance matrix used will be the Mahalanobis distance computed using the covariance matrix of the target distribution (weighted by \code{targetW} if supplied). If \code{"none"}, the distance matrix used will be the Euclidean distance matrix. Default is \code{"studentized"}. Abbreviations allowed.
}
  \item{lambda}{
a penalty on the sum of the squared weights. \code{lambda/nrow(sampleX)^2} times the sum of the squared weights is added to the objective function. Increasing \code{lambda} preserves the effective sample size of the weighted sample at the expense of balance.
}
}
\details{

}
\value{
A vector of weights for the treated group.
}

\author{
Noah Greifer
}
\seealso{
\code{\link{energybalance}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
