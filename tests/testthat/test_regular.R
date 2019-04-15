context("Testing rregular data sets")

test_that("regular", {
  sample_data <-
    c(0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0,
      1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1,
      1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1,
      0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0)
  M <- matrix(sample_data, ncol = 10)
  expect_equal(c(changes(M, 1:2)), c(0, 3, 3, 4, 10))
  expect_equal(c(changes(M, c(1, 3, 7))), c(7, 6, 3, 4, 10))
  expect_equal(sum(changes(M, 1:10)[1:4]), 90)
  expect_equivalent(round(rates(matrix(c(3, 4, 6, 7, 3), nrow = 1), CI = T)[2], digits = 2),
               1.71)
  expect_equivalent(rates(matrix(c(3, 4, 6, 7, 3), nrow = 1), CI = T)[7], 3)
  trans <- cetotrans(1, 1)
  expect_equal(round(lfun2(c = 1, e = 1, matrix(c(3, 4, 6, 7), nrow=1)),
                     digits = 2),
               round(log(trans[1]) * 7 + log(1 - trans[1]) * 13, digits = 2))
  expect_equal(nrow(regular_sampling_scheme(alonso15[[1]], 3:6, "Guild", 5,
                                            CI=T)), 7)
  expect_lt(abs(NLL_rss(alonso15[[1]], 3:6, 0.52, 0.39) -
                regular_sampling_scheme(alonso15[[1]], 3:6)$NLL), 0.1)

})
