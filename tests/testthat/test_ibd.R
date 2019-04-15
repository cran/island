### The following test doesn't works on CRAN, but works locally.

context("Testing IBD models")

# library(deSolve)
#
# test_that("IBD models conform to deterministic expectation", {
#   ### Preparing the function for deSolve
#   logistic <- function(time, state, pars){
#     r <- pars$r
#     K <- pars$K
#     mu <- pars$mu
#     N <- state
#     return(list(r * N * (1 - (N / K)) + mu))
#   }
#
#   # Parameters
#   r <- 0.1
#   mu <- 0.1
#   K <- 200
#
#   N <- 20
#
#   pars <- list(r = r, mu = mu, K = K)
#
#   set.seed(10110111)
#   # Solving
#   output <- ode(y = N, times = 0:150, func = logistic, method = "ode45", parms = pars)
#
#   df <- matrix(NA, 151, 50)
#   for (i in 1:50){
#     df[, i] <- ibd_models(20, 0.2, 0.1, 0.1, K = 200, time_v = 0:150,
#                           type = "Haegeman")[, 2]
#   }
#   # plot(0:150, rowMeans(df))
#   # lines(output)
#
#   expect_lt(output[151, 2] - rowMeans(df)[151], 0.05 * output[151, 2])
#
# })
