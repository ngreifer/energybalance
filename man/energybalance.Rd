\name{energybalance}
\alias{energybalance}

\title{Estimate Energy Balancing Weights}
\description{
Estimates energy balancing weights for matching a sample distribution to a target distribution.

}
\usage{
energybalance(sampleX, Z = NULL, targetX = NULL,
              sampleW = NULL, targetW = NULL,
              std = "studentized", improved = TRUE,
              lambda = 0)
}

\arguments{
  \item{sampleX}{
a matrix of covariates composing the sample distribution to be weighted.
}
  \item{Z}{
optional; a vector denoting treatment group membership. If non-\code{NULL}, each treatment group will be weighted to resemble the target distribution. If \code{NULL}, the entire sample distribution will be weighted to resemble the target distribution.
}
  \item{targetX}{
a matrix of covariates composing the target distribution. If \code{targetX} is \code{NULL} and \code{Z} is specified, \code{sampleX} will be used.
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
  \item{improved}{
\code{logical}; when \code{Z} is specified, whether to additional balance the distributions between treatment groups. If \code{TRUE}, the energy distance between treatment groups and between each treatment group and the target distribution will be minimized. If \code{FALSE}, only the energy distance between each treatment group and the target distribution will be minimized. Default is \code{TRUE} as recommended by Huling and Mak (2020). Ignored when \code{Z} is \code{NULL}.
}
  \item{lambda}{
a penalty on the sum of the squared weights. \code{lambda/nrow(sampleX)^2} times the sum of the squared weights is added to the objective function. Increasing \code{lambda} preserves the effective sample size of the weighted sample at the expense of balance.
}
}
\value{
A vector of weights, one for each row of \code{sampleX}.
}
\details{
\code{energybalance} is the main function of the \emph{energybalance} package. It estimates energy balancing weights for a broad set of scenarios, including estimating average treatment effects and generalizing or transporting estimates. The wrappers for \code{energybalance}--\code{\link{eb_ate}}, \code{\link{eb_target}}, and \code{\link{eb_mediation}}--simplify its use for specific purposes.

Essentially, \code{energybalance} estimates energy balancing weights to make the sample distribution (\code{sampleX}) resemble the target distribution (\code{targetX}) by minimizing the weighted energy distance between them. When a treatment vector (\code{Z}) is supplied, the weights make each treatment group resemble the target distribution. When sampling weights for the target distribution (\code{targetW}) are supplied, the weights make the sample resemble the weighted target distribution; this is especially useful when weighting a sample to resemble a representative sample that requires sampling weights.

The energy distance between two groups is dependent on the scale of the variables because its original formulation relies on the Euclidean distance matrix between the two groups. When \code{std} is \code{"studentized"} (the default), the scaled Euclidean distance is used instead, which eliminates the dependence on scale. When \code{std} is \code{"mahalanobis"}, the Mahalanobis distance is used instead, which eliminates the dependence on scale and the correlations between variables. Performance may vary for the different distances. Using the Mahalanobis distance eliminates the depence on the coding scheme used for factor variables (i.e., which category is dropped), but can yield poorer balancing performance for other covariates.

When \code{sampleW} are used, the final weights should be multiplied by \code{sampleW} prior to using them in an analysis.
}

\author{
Noah Greifer
}

\seealso{
The wrapper functions:
\itemize{
 \item {
  \code{\link{eb_ate}} for weights that estimate the average treatment effect (ATE).
 }
 \item {
  \code{\link{eb_target}} for weights that estimate the average treatment effect in the treated (ATT), allow generalizability or transportability of a sample to a target, or adjust for censoring by loss to follow-up.
 }
 \item {
  \code{\link{eb_mediation}} for weights that estimate the mean cross-world potential outcome for estimating natural direct and indirect effects in mediation analysis.
 }
}
}

\examples{

}
