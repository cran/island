% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataSimulation.R
\name{data_generation}
\alias{data_generation}
\alias{PA_simulation}
\title{Data simulation of colonization-extinction dynamics}
\usage{
data_generation(x, column, transitions, iter, times)

PA_simulation(x, column, transitions, times = 1)
}
\arguments{
\item{x}{A dataframe with the vector of initial absences and presences.}

\item{column}{A number indicating the column with the initial
presence-absence data.}

\item{transitions}{A matrix with the transition probabilities of the
simulation, in the form (T01, T10), that can contain one single pair or
multiple pairs.}

\item{iter}{Number of times that the specified dynamics should be repeated.}

\item{times}{Number of temporal steps to simulate.}
}
\value{
A matrix with species richness representing each row consecutive
  samples and each column a replica of the specified dynamics  or a matrix
  with presence-absence data for the specified dynamics, each row
  representing a species and each column consecutive samplings.
}
\description{
\code{data_generation} simulates species richness data according to the
stochastic model of island biogeography \cr \code{PA_simulation} simulates
presence-absence data according to the stochastic model of island
biogeography
}
\details{
To simulate community assembly, we need an initial vector of
presence-absence, from which the subsequent assembly process will be
simulated. This initial vector is considered as \code{x[, column]}.
}
\note{
You can simulate not only with a colonization and extinction pair, but
  with the pairs obtained from the environmental fit. In this case, you still
  have to indicate exactly the number of temporal steps that you are going to
  simulate.
}
\examples{
data_generation(as.data.frame(rep(0, 100)), 1,
matrix(c(0.5, 0.5), ncol = 2), 5, 25)
data_generation(alonso15[[1]], 3, matrix(c(0.5, 0.5), ncol = 2), 5, 25)
PA_simulation(as.data.frame(c(rep(0, 163), rep(1, 57))), 1, c(0.13, 0.19),
20)

}
\seealso{
\code{\link{cetotrans}} to obtain the transition probabilities
  associated with a colonization-extinction pair.
}
