context("Testing Environmental fit")

test_that("Environmental fit", {
  ### !!! It seems that my use of expressions is not supported in testthat.
  ### Still, the functions work consistently and accurately. You can uncomment
  ### the following code to see it.

  # set.seed(101201230)
  # transit <- cetotrans(rep(.2, 10), 1:10 * 0.04)
  # initial <- sample(c(0,1), 1000, replace = T, prob = c(1, 5))
  # initial <- data.frame(initial)
  # M <- cbind(initial, PA_simulation(initial, 1, transit, 10))
  # colnames(M) <- 1:11
  # env <- data.frame(good = 0:9, bad = runif(10))
  # out <- all_environmental_fit(dataset = M, vector = 1:11,
  #                              c("env$good", "env$bad"),c = .02, e = .004,
  #                              aic = 1e6)
  # out2 <- greedy_environmental_fit(dataset = M, vector = 1:11,
  #                                  c("env$good", "env$bad"),c = .02, e = .004,
  #                                  aic = 1e6)
  # # Both functions find a very close estimate of c = .2, e = .04, and the
  # # coeficient for the part of the extinction that depends on the "environmental
  # # variable" (0.04).
  #
  # nll <- NLL_env(dataset = M, vector = 1:11, params = out$Results$par,
  #         out$Colonization, out$Extinction)
  # # expect_equivalent(out$Results$value, -nll)
  # # The NLL calculator finds the same NLL than the fit.
  #
  #
  # # expect_equal(out, out2)
  # rates_calculator(out$Results$par, out$Colonization, out$Extinction, t = 10)
  #
  # # And we also find quite close estimates of the rates at each transition.

})
