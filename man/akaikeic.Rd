% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EvaluationCriteria.R
\name{akaikeic}
\alias{akaikeic}
\alias{akaikeicc}
\title{Akaike Information Criterion}
\usage{
akaikeic(NLL, k)

akaikeicc(NLL, k, n)
}
\arguments{
\item{NLL}{Negative Log-Likelihood of the model.}

\item{k}{Number of parameters of the model.}

\item{n}{Sample size.}
}
\value{
A number with the AIC value for a model with k parameters and
  negative log-likelihood NLL, or the AICc value for a model with k parameters,
  negative log-likelihood NLL and sample size n.
}
\description{
\code{akaikeic} calculates the Akaike Information Criterion (AIC) of a model.
\cr \code{akaikeicc} calculates the corrected Akaike Information Criterion
(AICc) for small samples.
}
\details{
\deqn{AIC = 2 * k + 2 * NLL} \deqn{AICc = 2 * k - 2 * lnL + 2 * k *
  (k + 1) / (n - k - 1)}
}
\examples{
akaikeic(1485.926, 3)
akaikeicc(736.47, 6, 15)
akaikeicc(736.47, 6, 100)
}
\seealso{
\code{\link{weight_of_evidence}}
}
