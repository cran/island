# This is the script for the functions and generics for the data generation
# tool.

#' Data simulation of colonization-extinction dynamics
#'
#' \code{data_generation} simulates species richness data according to the
#' stochastic model of island biogeography \cr \code{PA_simulation} simulates
#' presence-absence data according to the stochastic model of island
#' biogeography
#'
#' To simulate community assembly, we need an initial vector of
#' presence-absence, from which the subsequent assembly process will be
#' simulated. This initial vector is considered as \code{x[, column]}.
#'
#' @param x A dataframe with the vector of initial absences and presences.
#' @param column A number indicating the column with the initial
#'   presence-absence data.
#' @param transitions A matrix with the transition probabilities of the
#'   simulation, in the form (T01, T10), that can contain one single pair or
#'   multiple pairs.
#' @param iter Number of times that the specified dynamics should be repeated.
#' @param times Number of temporal steps to simulate.
#' @examples data_generation(as.data.frame(rep(0, 100)), 1,
#' matrix(c(0.5, 0.5), ncol = 2), 5, 25)
#' data_generation(alonso15[[1]], 3, matrix(c(0.5, 0.5), ncol = 2), 5, 25)
#' PA_simulation(as.data.frame(c(rep(0, 163), rep(1, 57))), 1, c(0.13, 0.19),
#' 20)
#'
#' @note You can simulate not only with a colonization and extinction pair, but
#'   with the pairs obtained from the environmental fit. In this case, you still
#'   have to indicate exactly the number of temporal steps that you are going to
#'   simulate.
#'
#' @return A matrix with species richness representing each row consecutive
#'   samples and each column a replica of the specified dynamics  or a matrix
#'   with presence-absence data for the specified dynamics, each row
#'   representing a species and each column consecutive samplings.
#'
#' @seealso \code{\link{cetotrans}} to obtain the transition probabilities
#'   associated with a colonization-extinction pair.
#'
#' @export
data_generation <- function(x, column, transitions, iter, times) {
  y <- data.frame()
  total <- matrix(NA, ncol = iter, nrow = times)

  if (nrow(transitions) >= 2) {
    T01 <- transitions[, 1]
    T10 <- transitions[, 2]
  } else {
    T01 <- rep(transitions[1], times)
    T10 <- rep(transitions[2], times)
  }

  for (i in 1:iter) {
    y <- cbind(x[, column])
    for (j in 1:times) {
      al <- stats::runif(nrow(y))
      for (k in 1:nrow(y)) {
        if (y[k, 1] == 0) {
          if (al[k] < T10[j])
            y[k, 1] <- 1
        }
        else {
          if (al[k] < T01[j])
            y[k, 1] <- 0
        }
      }
      total[j, i] <- colSums(y)
    }
  }
  total
}

#' @rdname data_generation
#' @export
PA_simulation <- function(x, column, transitions, times = 1) {

  y <- as.matrix(x[, column])


  if (is.null(nrow(transitions))) {
    T01 <- rep(transitions[1], times)
    T10 <- rep(transitions[2], times)

  } else {
    if(nrow(transitions) != 1){
      T01 <- c(transitions[,1])
      T10 <- c(transitions[,2])
      if(nrow(transitions) != times){
        warning("The number of transitions supplied is not the same as argument times. Using the number of transitions instead.",
                call. = F)
      }
      times <- nrow(transitions)
    } else {
      T01 <- rep(transitions[1], times)
      T10 <- rep(transitions[2], times)
    }
  }
  result <- matrix(NA, ncol = times, nrow = nrow(x))
  result[, 1] <- y

  for (j in 1:times) {
    al <- stats::runif(nrow(x))
    for (k in 1:nrow(x)) {
      if (y[k, 1] == 0) {
        if (al[k] < T10[j]) {
          y[k, 1] <- 1
        }
      }
      else {
        if (al[k] < T01[j]) {
          y[k, 1] <- 0
        }
      }
    }
    result[, j] <- y
  }
  result
}
