\name{eb_target}
\alias{eb_target}
\title{
Estimate Energy Balancing Weights To Match a Sample To a Target
}
\description{
\code{eb_target} is a wrapper for \code{\link{energybalance}} that focuses on using it to estimate energy balancing weights for matching a sample distribution to a target distribution, i.e., for the purposes of generalizability, transportability, or estimating the average treatment effect in the treated (ATT). The energy distance between the sample distribution and the target distribution is minimized.
}
\usage{
eb_target(sampleX, targetX, sampleW = NULL,
          targetW = NULL, std = "studentized")
}
\arguments{
  \item{sampleX}{
a matrix of covariates composing the sample distribution to be weighted.
}
  \item{targetX}{
a matrix of covariates composing the target distribution.
}
  \item{sampleW}{
optional; a vector of sampling weights for the sample. The product of the estimated weights and \code{sampleW} will be the distribution-matching weights.
}
  \item{targetW}{
optional; a vector of sampling weights for the target.
}
  \item{std}{
\code{character}; whether to standardize the covariates. If \code{"studentized"}, the distance matrix used will be the Euclidean distance with each variable scaled using its standard deviation in the target distribution (weighted by \code{targetW} if supplied). If \code{"mahalanobis"}, the distance matrix used will be the Mahalanobis distance computed using the covariance matrix of the target distribution (weighted by \code{targetW} if supplied). If \code{"none"}, the distance matrix used will be the Euclidean distance matrix. Default is \code{"studentized"}. Abbreviations allowed.
}
}
\value{
A vector of weights, one for each row of \code{sampleX}.
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
