context("Testing irregular data sets")

test_that("single dataset", {
  sample_data <-
  c(0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0,
    1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1,
    1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1,
    0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0)
  M <- matrix(sample_data, ncol = 10)
  colnames(M) <- 1:10
  expect_equivalent(c(colSums(incounts(M, 1:10))), c(changes(M, 1:10)[1:4]))
  expect_equivalent(c(colSums(incounts(M, c(3, 5, 9)))),
                    c(changes(M, c(3, 5, 9))[1:4]))
  expect_equivalent(c(colSums(incounts(M, 1:10, assembly = T))),
              c(changes(M, 1:10)[1:4]) + c(0, sum(M[, 1]),
                                           nrow(M) - sum(M[, 1]), 0))

  expect_equal(round(irregular_single_dataset(M, 1:10, c = 1, e = 1, CI = T),
                     1), round(regular_sampling_scheme(M, 1:10, CI = T), 1))

  expect_equal(round(irregular_single_dataset(M, 1:10, c = 1, e = 1, jacobian = T,
                                              CI = T), 1),
               round(regular_sampling_scheme(M, 1:10, CI = T), 1))

  expect_equal_to_reference(irregular_single_dataset(simberloff[[2]], 3:17,
                          column = "Tax. Unit 1", c = 0.00000001, e = 0.0000001,
                          n = 15, assembly = T, jacobian = T, CI = T), "i1.rds")
  expect_lt(sum(irregular_single_dataset(simberloff[[2]], 3:17,
                column = "Tax. Unit 1", c = 0.00000001, e = 0.0000001, n = 15,
                assembly = T, jacobian = T, CI = T)[, 4] -
                irregular_single_dataset(simberloff[[2]], 3:17,
                column = "Tax. Unit 1", c = 0.00000001, e = 0.0000001, n = 15,
                step = 0.0001, assembly = T, jacobian = T, CI = T)[, 4]), 0.001)
  expect_equal(colnames(irregular_single_dataset(simberloff[[2]], 3:17,
                        column = "Tax. Unit 1", c = 0.00000001, e = 0.0000001,
                        n = 15, assembly = T, jacobian = T, CI = T)),
               colnames(irregular_single_dataset(simberloff[[2]], 3:17,
                        column = "Tax. Unit 1", c = 0.00000001, e = 0.0000001,
                        n = 15, step = 0.001, assembly = T, jacobian = T, CI = T)))
  expect_lt(abs(irregular_single_dataset(simberloff[[1]], 3:17, 0.001, 0.001)$NLL -
                  NLL_isd(simberloff[[1]], 3:17, 0.0038, 0.0086)), 0.1)
})

test_that("multiple datasets", {
  expect_equal(sum(irregular_multiple_datasets(list(simberloff[[1]]), list(3:17), 0.00000001, 0.00000001,
            verbose = F, jacobian = T, assembly=T, CI = T) -
  irregular_single_dataset(simberloff[[1]], 3:17, c = 0.00000001, e = 0.0000001,
                           assembly = T, jacobian = T, CI = T) > 1e-6), 0)
  expect_equal_to_reference(wrapper2(simberloff, list(3:17, 3:18, 3:17, 3:19, 3:17, 3:16), 0.001, 0.001,
            verbose = F, jacobian = T, assembly=T, CI = T, step = NULL), "i2.rds")
  expect_equal(sum((irregular_multiple_datasets(list(simberloff[[2]]), list(3:17), 0.001, 0.001, "Tax. Unit 1", n = 25, CI = TRUE)[-1] -
  irregular_single_dataset(simberloff[[2]], 3:17,
                           column = "Tax. Unit 1", c = 0.00000001, e = 0.0000001, n = 25,
                           CI = T)[-1]) > 2e-5), 0)
  expect_lt(abs(irregular_multiple_datasets(simberloff, list(3:17, 3:18, 3:17,
            3:19, 3:17, 3:16), 0.001, 0.001)$NLL - NLL_imd(simberloff, list(3:17,
            3:18, 3:17, 3:19, 3:17, 3:16), 0.0051, 0.0117)), 0.1)
  })
