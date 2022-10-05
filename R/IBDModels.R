#' Inmigration, birth, death- models
#'
#' \code{ibd_models} simulates population dynamics under three different
#' inmigration, birth and death models.
#'
#' We have included three different stochastic models: Kendall (1948) seminal
#' work, Alonso & McKane (2002) mainland-island model, and Haegeman & Loreau
#' (2010) basic population model with denso-dependent deaths. These models are
#' different formulations of a population dynamics with three basic processes:
#' birth, death and inmigration of individuals. For the specifics, please refer
#' to the original articles.
#'
#' @param n0 Initial number of individuals in the population.
#' @param beta Birth rate, in time^(-1) units.
#' @param delta Death rate, in time^(-1) units.
#' @param mu Inmigration rate, in time^(-1) units.
#' @param K Carrying capacity.
#' @param time_v Vector of times to sample. Must start with 0.
#' @param type Type of inmigration, birth, death- model used to simulate the
#'   dynamics. This must be \code{"Kendall"}, \code{"Alonso"} or \code{"Haegeman"}.
#' @examples ibd_models(n0 = 0, beta = 0.4, delta = 0.3, mu = 0.2,
#' time_v = 0:20, type = "Kendall")
#' ibd_models(n0 = 0, beta = 0.4, delta = 0.3, mu = 0.1, K = 30,
#' time_v = 0:20, type = "Alonso")
#'
#' @note Haegeman & Loreau model specification breaks for high values of
#'   \code{n0} when the birth rate is lower than the death rate.
#'
#' @return A data.frame with two columns: one with the time vector and the other
#'   with the number of individuals at those times.
#'
#' @references Kendall, D. G. (1948). On some modes of population growth leading
#'   to R. A. Fishers logarithmic series distribution. \emph{Biometrika},
#'   \bold{35}, 6--15. \cr \cr Haegeman, B. and Loreau, M. (2010). A
#'   mathematical synthesis of niche and neutral theories in community ecology.
#'   \emph{Journal of Theoretical Biology}, \bold{269(1)}, 150--165. \cr \cr
#'   Alonso, D. and McKane, A (2002). Extinction Dynamics in Mainland--Island
#'   Metapopulations: An N -patch Stochastic Model. \emph{Bulletin of
#'   Mathematical Biology}, \bold{64}, 913--958.
#'
#' @export
ibd_models <- function(n0, beta, delta, mu, K = NULL, time_v, type){
  # Select type of model
  if (type == "Kendall") {
    gn <- function(n) { beta * n + mu }
    rn <- function(n) { delta * n }
  } else {
    if (type == "Alonso") {
      gn <- function(n) { beta * n * (1 - (n / K)) + mu * (K - n) }
      rn <- function(n) { delta * n }
    } else {
      if (type == "Haegeman") {
        gn <- function(n) { beta * n + mu }
        rn <- function(n) { delta * n  + ((beta - delta) / K) * (n ^ 2)}
        if (delta > beta){
          nmax <- K * (- delta / (beta - delta))
          cat("If n goes above ", nmax, ", the function will stop.")
        }
      } else {
        stop("Type of model incorrectly specified.")
      }
    }
    if (is.null(K)) stop("Carrying capacity K not specified.")
  }

  # Initialize variables
  n <- n0
  t <- 0
  out <- matrix(ncol = 2, nrow = length(time_v))
  out[1,] <- c(0, n0)

  # Start simulations
  for(i in 2:length(time_v)){
    while (t < time_v[i]){
      n0 <- n

      # Rates
      r.gn <- gn(n)
      r.rn <- rn(n)

      # Sampling
      dt <- stats::rexp(1, r.gn + r.rn)
      event <- sample(c(1, -1), 1, prob = c(r.gn, r.rn))
      n <- n + event
      t <- t + dt
    }
    out[i,] <- c(time_v[i], n0)
  }
  colnames(out) <- c("Time", "n")
  out
}

