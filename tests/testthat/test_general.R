context("Testing General")

test_that("cetotrans", {
  expect_equal(nrow(cetotrans(rep(0.1, 3), rep(0.2, 3))), 3)
  expect_equal_to_reference(cetotrans(0.1, 0.2), "g1.rds")
  expect_equal_to_reference(cetotrans(rep(0.1, 3), rep(0.2, 3), 1:3), "g2.rds")
})
test_that("Null Models", {
  set.seed(10110111)
  idaho.sim <- data_generation(as.data.frame(c(rep(0, 163), rep(1, 57))),
                               1,
                               matrix(c(0.162599, 0.111252), ncol = 2),
                               250,
                               20)
  idaho.me <- c(57, apply(idaho.sim, 1, quantile, 0.5))
  t1 <- r_squared(colSums(idaho[[1]][,3:23]), idaho.me, 220)
  t2 <- null_model(colSums(idaho[[1]][,3:23]), 220)
  t3 <- simulated_model(colSums(idaho[[1]][,3:23]), idaho.me)

  expect_equal_to_reference(list(t1 = r_squared(colSums(idaho[[1]][,3:23]), idaho.me, 220),
                                 t2 = null_model(colSums(idaho[[1]][,3:23]), 220),
                                 t3 = simulated_model(colSums(idaho[[1]][,3:23]), idaho.me)),
                            "g3.rds")

})

test_that("Akaike", {
  expect_equal(akaikeic(710.79, 42), 1505.58)
  expect_equal(akaikeicc(710.79, 42, 1248), 1508.578, tol = 0.001)

  models <- c("Best_3k", "Best_4k", "Best_5k", "Best_6k", "Best_7k",
    "Best_8k", "Best_9k")

  aks <- c(2977.852, 2968.568, 2957.384, 2952.618,
    2949.128, 2947.038, 2943.480)

  expect_equal_to_reference(weight_of_evidence(cbind(models, aks)), "g4.rds")
})

test_that("Simulations", {
  set.seed(10110111)
  aaa <- data_generation(as.data.frame(rep(0, 100)),
                         1,
                         matrix(c(0.5, 0.5), ncol = 2),
                         5,
                         25)
  expect_equal(ncol(aaa), 5)
  expect_equal(nrow(aaa), 25)

  set.seed(10110111)
  bbb <- colSums(PA_simulation(as.data.frame(rep(0, 100)),
                1,
                matrix(c(0.5, 0.5), ncol = 2),
                times = 25))

  expect_equal(aaa[, 1], bbb)
})
