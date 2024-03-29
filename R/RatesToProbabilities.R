#' From rates to probabilities
#'
#' \code{cetotrans} calculates transition probabilities from colonization and
#' extinction rates for a determined interval of time, when provided.
#'
#' @param c Colonization rate.
#' @param e Extinction rate.
#' @param dt Interval of time or a vector of time intervals.
#'
#' @details Given a pair of colonization and extinction rates, we can calculate
#'   the transition probabilities with the following equations: \deqn{T_{01} = (e
#'   / (c + e)) * (1 - exp( - (c + e) * dt))} \deqn{T_{10} = (c / (c + e)) * (1 -
#'   exp( - (c + e) * dt))}
#'
#' @return A matrix with the transition probabilities \eqn{T_{01}} and \eqn{T_{10}} of the Markov chain
#'   associated with the specified colonization and extinction rates.
#' @examples cetotrans(0.13, 0.19)
#' cetotrans(0.2, 0.2, 2)
#' @export

cetotrans <- function(c, e, dt = 1){

  ### This function calculates T01 and T10 given a c, e, and dt.

  T01 <- (e / (c + e)) * (1 - exp( - (c + e) * dt))
  T10 <- (c / (c + e)) * (1 - exp( - (c + e) * dt))
  out <- matrix(c(T01, T10), ncol = 2)
  colnames(out) <- c("T01", "T10")
  out
}
