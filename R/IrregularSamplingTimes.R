#' c/e rates for irregular samplings in a dataset
#'
#' \code{irregular_single_dataset} estimates colonization and extinction rates
#' in a single dataset with irregular sampling scheme. \cr \code{NLL_isd} returns
#' the Negative Log-Likelihood of a pair of colonization and extinction rates
#' for an irregular sampling scheme in a single dataset.
#'
#' @param dataframe A single dataframe.
#' @param vector A vector indicating the columns with presence-absence data.
#' @param c Tentative colonization rate.
#' @param e Tentative extinction rate.
#' @param column The name of the column with groups to calculate their c_e pair.
#' @param n Minimal number of rows for each group
#' @param step Accuracy to calculate the c_e pairs with.
#' @param assembly Logical indicating if the assembly starts from zero species
#'   or not.
#' @param jacobian Logical. Use the semianalytical method to estimate colonization
#'   and extinction rates?
#' @param verbose Logical. If TRUE, gives the output of the optimizer or the
#'   numerical solver that finds the values of c and e.
#' @param CI Logical. If TRUE, gives the confidence interval of the colonization
#'   and extinction rates.
#' @export
#' @seealso \code{\link{regular_sampling_scheme}},
#'   \code{\link{irregular_multiple_datasets}}
#' @note The columns with the presence-absence data should have the day of that
#'   sampling on the name of the column in order to calculate colonization and
#'   extinction.
#' @examples irregular_single_dataset(simberloff[[1]], 3:17, 0.001, 0.001)
#' irregular_single_dataset(simberloff[[1]], 3:17, column = "Tax. Unit 1",
#' 0.001, 0.001, 3)
#' \dontrun{
#' irregular_single_dataset(simberloff[[1]], 3:17, column = "Tax. Unit 1",
#' 0.001, 0.001, 3, 0.00001)
#' }
#' NLL_isd(simberloff[[1]], 3:17, 0.0038, 0.0086)
#' @return \code{irregular_single_dataset} returns a dataframe with colonization
#'   and extinction rates and their upper and lower confidence interval, and if
#'   needed, the names of the groups to which colonization and extinction rates
#'   have been calculated. \code{NLL_isd} gives the NLL for a single dataset in
#'   an irregular sampling scheme given a specific c and e.
irregular_single_dataset <- function(dataframe, vector, c, e, column = NULL,
            n = NULL, step = NULL, assembly = F, jacobian = F, verbose = F, CI=F) {

  if (is.null(column)) {
    wrapper(dataset = dataframe, vector = vector, c = c, e = e, step = step,
            assembly = assembly, jacobian = jacobian, verbose = verbose, CI = CI)
  } else {
    levelwrap2(dataframe, vector, column, c, e, n, step, assembly, jacobian,
               verbose, CI)
  }
}

#' @rdname irregular_single_dataset
#' @export

NLL_isd <- function(dataframe, vector, c, e, assembly = F){

  co <- incounts(dataframe, vector, assembly)
  t <- times(dataframe, vector, assembly)

  -dtlsolver(c(c, e), t = t, n00 = co[, 3], n10 = co[, 2], n01 = co[, 1],
                    n11 = co[, 4])
}


#' c/e rates for irregular samplings in multiple datasets
#'
#' \code{irregular_multiple_datasets} estimates colonization and extinction
#' rates for data in several datasets. \cr \code{NLL_imd} returns the Negative
#' Log-Likelihood of a pair of colonization and extinction rates for irregular
#' sampling schemes in several single dataset.
#' @param list A list of dataframes.
#' @param vectorlist A list of vectors indicating the columns with
#'   presence-absence data.
#' @param c Tentative colonization rate.
#' @param e Tentative extinction rate.
#' @param column The name of the column with groups to calculate their c_e pair.
#' @param n Minimal number of rows for each group.
#' @param step Accuracy to calculate the c_e pairs with.
#' @param assembly Logical indicating if the assembly starts from zero species
#'   or not.
#' @param jacobian Logical. Use the semianalytical method to estimate colonization
#'   and extinction rates?
#' @param verbose Logical. If TRUE, gives the output of the optimizer or the
#'   numerical solver that finds the values of c and e.
#' @param CI Logical. If TRUE, gives the confidence interval of the colonization
#'   and extinction rates.
#'
#' @export
#' @seealso \code{\link{regular_sampling_scheme}},
#'   \code{\link{irregular_single_dataset}}
#' @note The columns with the presence-absence data should have the day of that
#'   sampling on the name of the column in order to calculate colonization and
#'   extinction.
#' @examples irregular_multiple_datasets(simberloff, list(3:17, 3:18, 3:17,
#' 3:19, 3:17, 3:16), 0.001, 0.001)
#' \dontrun{
#' irregular_multiple_datasets(simberloff, list(3:17, 3:18, 3:17, 3:19, 3:17,
#'  3:16), 0.001, 0.001, "Tax. Unit 1", n = 13)
#' irregular_multiple_datasets(simberloff, list(3:17, 3:18, 3:17, 3:19, 3:17,
#'  3:16), 0.001, 0.001, "Tax. Unit 1", n = 13, CI = TRUE)
#'  }
#'  NLL_imd(simberloff, list(3:17, 3:18, 3:17, 3:19, 3:17, 3:16), 0.0051, 0.0117)
#' @return \code{irregular_multiple_datasets} returns a dataframe with
#'   colonization and extinction rates and their upper and lower confidence
#'   interval, and if needed, the names of the groups to which colonization and
#'   extinction rates have been calculated. \code{NLL_imd} gives the NLL for a
#'   multiple datasets with irregular sampling schemes given a specific c and e.

irregular_multiple_datasets <- function(list, vectorlist, c, e, column = NULL,
           n = NULL, step = NULL, assembly = F, jacobian = F, verbose = F, CI = F) {

  if (is.null(column)) {
    wrapper2(list, vectorlist, c, e, step, assembly, jacobian, verbose, CI)
  } else {
    levelwrapb(list, vectorlist, column, c, e, n, step, assembly, jacobian, verbose,
               CI)
  }
}

#' @rdname irregular_multiple_datasets
#' @export

NLL_imd <- function(list, vectorlist, c, e, assembly = F){

  time <- vector("list", length(list))
  abs  <- vector("list", length(list))
  pre  <- vector("list", length(list))
  col  <- vector("list", length(list))
  ext  <- vector("list", length(list))

  for (i in 1:length(list)) {
    dataset <- list[[i]]
    vector <- vectorlist[[i]]
    co <- data.frame()
    co <- incounts(dataset, vector, assembly)
    time[i] <- data.frame(times(dataset, vector, assembly))
    abs[i]  <- data.frame(co[, 3])
    pre[i]  <- data.frame(co[, 4])
    col[i]  <- data.frame(co[, 2])
    ext[i]  <- data.frame(co[, 1])

  }
    - dtlsolver2(c(c, e), t = time, abs = abs, col = col, ext = ext,
                     pre = pre)
}

dtlsolver <- function(x, t, n00, n10, n01, n11) {

  ### This function calculates the sum of the loglikelihood for a set x
  ### with two components, c and e, given a time series t and the number of
  ### transitions n00, n10, n01 and n11 in each time. Necessary for
  ### the solver.

  c   <- x[1]
  e   <- x[2]
  sum <- 0

  for (i in 1:length(t)) {
    ll <- n00[i] * log(1 - (c / (e + c)) * (1 - exp( - (e + c) * t[i]))) +
          n11[i] * log(1 - (e / (e + c)) * (1 - exp( - (e + c) * t[i]))) +
          n10[i] * log( (c / (e + c)) * (1 - exp( - (e + c) * t[i]))) +
          n01[i] * log( (e / (e + c)) * (1 - exp( - (e + c) * t[i])))
    sum <- sum + ll
  }
  sum
}

dtlsolver2 <- function(x, t, abs, col, ext, pre) {

  ### Another version of the dtlsolver. Posibly, made to work with multiple
  ### discontinous times.

  c   <- x[1]
  e   <- x[2]
  res <- 0

  for (i in 1:length(t)) {
    if(is.null(t[[i]])) next
    sum   <- 0
    times <- t[[i]]
    n00   <- abs[[i]]
    n10   <- col[[i]]
    n01   <- ext[[i]]
    n11   <- pre[[i]]
    for (j in 1:length(times)) {
      ll <- n00[j] * log(1 - (c / (e + c)) * (1 - exp( - (e + c) * times[j]))) +
            n11[j] * log(1 - (e / (e + c)) * (1 - exp( - (e + c) * times[j]))) +
            n01[j] * log( (e / (e + c)) * (1 - exp( - (e + c) * times[j]))) +
            n10[j] * log( (c / (e + c)) * (1 - exp( - (e + c) * times[j])))
      sum <- sum + ll
    }
    res <- res + sum
  }
  res
}

incounts <- function(dataset, vector, assembly = FALSE) {
  x <- dataset
  resultado <- data.frame()
  ### This function returns, from a dataset with data of presence-absence
  ### in each column of vector and the next in the dataset, a dataframe
  ### with four columns, corresponding to N01, N10, N00 and N11 respectively,
  ### and each row corresponding with a different temporal change. This
  ### function works for the cases in which the first column of the vector
  ### matches with the presence-absence at time 0 and those times when the
  ### first column is not at time 0, if the names of the columns are numeric.
  if (assembly) {

    if (as.numeric(colnames(x)[vector])[1] != 0) {
      resultado <- c(0, sum(dataset[, vector[1]]),
                     (nrow(x) - sum(dataset[, vector[1]])), 0)
    }

    for (j in 1:(length(vector) - 1)) {
      N10 <- 0; N01 <- 0; N00 <- 0; N11 <- 0
      for (i in 1:nrow(x)) {
        if (x[i, vector[j]] < x[i, vector[j + 1]]) N10 <- N10 + 1
        if (x[i, vector[j]] > x[i, vector[j + 1]]) N01 <- N01 + 1
        if (x[i, vector[j]] == x[i, vector[j + 1]] && x[i, vector[j]] == 0) {
          N00 <- N00 + 1
        }
        if (x[i, vector[j]] == x[i, vector[j + 1]] && x[i, vector[j]] == 1) {
          N11 <- N11 + 1
          }
      }
      resultado <- rbind(resultado, c(N01, N10, N00, N11))
    }
    rownames(resultado) <- NULL
    resultado
  } else {
    for (j in 1:(length(vector) - 1)){
      N10 <- 0; N01 <- 0; N00 <- 0; N11 <- 0
      for (i in 1:nrow(x)){
        if (x[i, vector[j]] < x[i, vector[j + 1]]) N10 <- N10 + 1
        if (x[i, vector[j]] > x[i, vector[j + 1]]) N01 <- N01 + 1
        if (x[i, vector[j]] == x[i, vector[j + 1]] && x[i, vector[j]] == 0) {
          N00 <- N00 + 1
          }
        if (x[i, vector[j]] == x[i, vector[j + 1]] && x[i, vector[j]] == 1) {
          N11 <- N11 + 1
          }
      }
      resultado <- rbind(resultado, c(N01, N10, N00, N11))
    }
    rownames(resultado) <- NULL
    resultado
  }
}

index <- function(list, column) {
  out <- vector()
  for (i in 1:length(list)) {
    x <- list[[i]]
    index <- levels (factor(x[, column]))
    out <- c(out, index)
  }
  out <- unique(levels(factor(out)))
  out
}

levelwrap <- function(list, vectorlist, column, c, e, n, assembly) {

  ### This function wraps the process of obtaining c and e in multiple
  ### datasets with different times (with different intervals) for the
  ### different levels contained in the specified "column". It also gives the
  ### C.I. calculated with the hessian.

  groups <- index(list, column)
  out2 <- data.frame()

  for (k in 1:length(groups)) {
    min  <- 0
    time <- vector("list", length(list))
    abs  <- vector("list", length(list))
    pre  <- vector("list", length(list))
    col  <- vector("list", length(list))
    ext  <- vector("list", length(list))

    for (i in 1:length(list)) {
      x <- list[[i]]
      dataset <- x[x[, column] == groups[k], ]
      if (nrow(dataset) == 0) next
      min <- min + nrow(dataset)
      vector <- vectorlist[[i]]
      co <- data.frame()
      co <- incounts(dataset, vector, assembly)
      time[i] <- data.frame(times(dataset, vector, assembly))
      abs[i]  <- data.frame(co[, 3])
      pre[i]  <- data.frame(co[, 4])
      col[i]  <- data.frame(co[, 2])
      ext[i]  <- data.frame(co[, 1])
    }

    if (min < n) next
    r <- stats::optim(c(c, e), dtlsolver2, t = time, abs = abs, col = col,
         ext = ext, pre = pre, control = list(maxit = 5000000, fnscale = - 1),
         hessian = T)

    fisher_info <- solve(-r$hessian)

    prop_sigma  <- sqrt(diag(fisher_info))
    prop_sigma  <- diag(prop_sigma)
    upper <- c(r$par[1] + 1.96 * prop_sigma[1, 1],
               r$par[2] + 1.96 * prop_sigma[2, 2])
    lower <- c(r$par[1] - 1.96 * prop_sigma[1, 1],
               r$par[2] - 1.96 * prop_sigma[2, 2])
    interval <- data.frame(value = r$par, upper = upper, lower = lower)
    out <- data.frame(groups[k], interval[1, ], interval[2, ])
    colnames(out) <- c("Group", "c", "Up-c", "Lo-c", "e", "Up-e", "Lo-e")
    out2 <- rbind(out2, out)
  }
  out2
}

levelwrap2 <- function(dataset, vector, column, c, e, n, step = NULL, assembly,
                       jacobian, verbose, CI) {

  ### This function wraps the process of obtaining c and e in a dataset with
  ### different times (with different intervals) for the different levels
  ### contained in the specified "column". It also gives the C.I. calculated
  ### with the hessian.

  groups <- levels(factor(dataset[, column]))
  out2   <- data.frame()
  out    <- data.frame()

  for (k in 1:length(groups)) {
    x  <- dataset[dataset[, column] == groups[k], ]
    if(nrow(x) < n) next

    outb <- wrapper(x, vector, c, e, step, assembly, jacobian, verbose, CI,
                    group = groups[k])

    out <- data.frame(groups[k], outb)
    colnames(out) <- c("Group", colnames(outb))
    out2 <- rbind(out2, out)
  }
  out2
}

levelwrapb <- function(list, vectorlist, column, c, e, n, step = NULL, assembly,
                       jacobian, verbose, CI) {

  ### This function wraps the process of obtaining c and e in multiple
  ### datasets with different times (with different intervals) for the
  ### different levels contained in the specified "column".

  groups <- index(list, column)
  out2   <- data.frame()

  for (k in 1:length(groups)) {
    min  <- 0
    list2 <- vector("list", length(list))
    for (i in 1:length(list)) {
      x <- list[[i]]
      dataset <- x[x[, column] == groups[k], ]
      if (nrow(dataset) == 0) next
      min <- min + nrow(dataset)
      list2[[i]] <- dataset
    }

    if (min < n) next
    empty <- !sapply(list2, is.null)
    list2 <- list2[empty]
    vectorlist2 <- vectorlist[empty]

    outb <- wrapper2(list2, vectorlist2, c, e, step, assembly, jacobian, verbose,
                     CI, group = groups[k])
    outb[7] <- min

    out <- data.frame(groups[k], outb)

    colnames(out) <- c("Group", colnames(outb))
    out2 <- rbind(out2, out)
  }
  out2
}

times <- function(dataset, vector, assembly = FALSE) {

  ### Given a dataset, it calculates the changes in time from column to
  ### column in the vector. The names of the columns need to be numbers. The
  ### output is a vector of increments in time from column to column.
  if (assembly) {
    times <- as.numeric(colnames(dataset)[vector])
    if (times[1] != 0) {
      out <- c(1:length(times))
    } else {
        out <- c(1:(length(times) - 1))
        }
    if (times[1] != 0) {
      out[1] <- times[1]
    } else {
        out[1] <- times[2] - times[1]
        }
    for (i in 2:length(out)) {
      if (times[1] == 0) {
        out[i] <- times[i + 1] - times[i]
      } else {
        out[i] <- times[i] - times[i - 1]
      }
    }
    out
  } else {
    times <- as.numeric(colnames(dataset)[vector])
    out <- c(1:(length(times) - 1))
    for (i in 1:length(out)) {

      out[i] <- times[i + 1] - times[i]
    }
    out
  }
}

wrapper <- function(dataset, vector, c, e, step, assembly, jacobian, verbose, CI,
                    group = "not defined") {

  ### This function wraps the whole process of obtaining c and e from a
  ### dataset with discontinuous times.

  co <- incounts(dataset, vector, assembly)
  t <- times(dataset, vector, assembly)

  if(jacobian){
    if (!requireNamespace("rootSolve", quietly = TRUE)) {
      stop("Package rootSolve needed for this function to work. Please install it.",
           call. = FALSE)
    }

    r <- rootSolve::multiroot(unequal_jacobian, c(c, e), N00 = co[,3],
                              N01 = co[,1], N10 = co[,2], N11 = co[,4], dt = t,
                              positive = T)

    if(is.nan(r$estim.precis)){
      stop(paste("Numerical solver didn't converge for", group,
                 " group. Please use lower priors for c and e or exclude it from the analysis."),
           call. = F)
    }
    if(r$estim.precis > 0.0001){
      warning(paste0("Numerical solver probably converged to a local optimum in group ",
                     group,
                     ". Please use lower priors for c and e or exclude it from the analysis."),
              call. = F)
    }

    c <- r$root[1]
    e <- r$root[2]
    ll <- dtlsolver(c(c, e), t = t, n00 = co[, 3], n10 = co[, 2], n01 = co[, 1],
                    n11 = co[, 4])
  } else {
    r <- stats::optim(c(c, e), dtlsolver, t = t, n00 = co[, 3], n10 = co[, 2],
                      n01 = co[, 1], n11 = co[, 4],
                      control = list(maxit = 5000000, fnscale = - 1),
                      hessian = F)

    c <- r$par[1]
    e <- r$par[2]
    ll <- r$value
  }
  if(verbose) print(r)
  if(CI){
    if(is.null(step)){
      while(CI){
        mhh <- -matrix(c(sum(hessian_c_r(n00 = co[3], n11 = co[4], n01 = co[1],
                                         n10 = co[2], dt = t, c = c, e = e)),
                         0,
                         0,
                         sum(hessian_e_r(n00 = co[3], n11 = co[4], n01 = co[1],
                                         n10 = co[2], dt = t, c = c, e = e))),
                       nrow = 2)

        h <- suppressWarnings(try(sqrt(diag(solve(mhh))) * 1.96, silent = T))
        if (inherits(h, "try-error")) {
          clo <- NA
          cup <- NA
          elo <- NA
          eup <- NA
          break
        } else {

          if (is.nan(h[1])) h[1] <- .99 * c
          if (is.nan(h[2])) h[2] <- .99 * e

          lower <- c(c - h[1], e - h[2])
          upper <- c(c + h[1], e + h[2])
          n <- 0
          if (lower[1] < 0) {
            fff <- c
            dx <- c
          } else {
            fff <- lower[1]
            f.fff <- dtlsolver(c(fff, e), t = t, n00 = co[, 3], n10 = co[, 2],
                               n01 = co[, 1], n11 = co[, 4])
            dx <- fff - c
            lld1.96 <- ll - f.fff < 1.96
            while (lld1.96) {
              fff <- fff + dx
              if (fff < 0) {
                fff <- fff - dx
                dx <- fff / 2
                fff <- fff - dx
              }
              f.fff <- dtlsolver(c(fff, e), t = t, n00 = co[, 3], n10 = co[, 2],
                                 n01 = co[, 1], n11 = co[, 4])
              n <- n + 1
              if (n > 1000) {
                break
              }
              lld1.96 <- ll - f.fff < 1.96
            }
          }
          dx <- -dx
          tol <- F
          while (!tol) {
            dx <- dx / 2
            fff <- fff + dx
            if (fff < 0) {
              fff <- fff - dx
              dx <- fff / 2
              fff <- fff - dx
            }
            f.fff <- dtlsolver(c(fff, e), t = t, n00 = co[, 3], n10 = co[, 2],
                               n01 = co[, 1], n11 = co[, 4])
            lld <- ll - f.fff
            if (lld < 1.96) dx <- -abs(dx) else dx <- abs(dx)
            n <- n + 1
            if (n > 1000) {
              fff <- NA
              break
            }
            tol <- abs(1.96 - lld) < 1e-3

          }
          clo <- fff

          n <- 0
          if (lower[2] < 0) {
            fff <- e
            dx <- e
          } else {
            fff <- lower[2]
            f.fff <- dtlsolver(c(c, fff), t = t, n00 = co[, 3], n10 = co[, 2],
                               n01 = co[, 1], n11 = co[, 4])
            dx <- fff - e
            lld1.96 <- ll - f.fff < 1.96
            while (lld1.96) {
              fff <- fff + dx
              if (fff < 0) {
                fff <- fff - dx
                dx <- fff / 2
                fff <- fff - dx
              }
              f.fff <- dtlsolver(c(c, fff), t = t, n00 = co[, 3], n10 = co[, 2],
                                 n01 = co[, 1], n11 = co[, 4])
              n <- n + 1
              if(n > 1000){
                break
              }
              lld1.96 <- ll - f.fff < 1.96
            }
          }
          dx <- -dx
          tol <- F
          while (!tol) {
            dx <- dx / 2
            fff <- fff + dx
            if (fff < 0) {
              fff <- fff - dx
              dx <- fff / 2
              fff <- fff - dx
            }
            f.fff <- dtlsolver(c(c, fff), t = t, n00 = co[, 3], n10 = co[, 2],
                               n01 = co[, 1], n11 = co[, 4])
            lld <- ll - f.fff
            if (lld < 1.96) dx <- -abs(dx) else dx <- abs(dx)
            if (n > 1000) {
              fff <- NA
              break
            }
            tol <- abs(1.96 - lld) < 1e-3
          }
          elo <- fff

          n <- 0
          fff <- upper[1]
          f.fff <- dtlsolver(c(fff, e), t = t, n00 = co[, 3], n10 = co[, 2],
                             n01 = co[, 1], n11 = co[, 4])
          dx <- fff - c
          lld1.96 <- ll - f.fff < 1.96
          while (lld1.96) {
            fff <- fff + dx
            f.fff <- dtlsolver(c(fff, e), t = t, n00 = co[, 3], n10 = co[, 2],
                               n01 = co[, 1], n11 = co[, 4])
            n <- n + 1
            if (n > 1000) {
              break
            }
            lld1.96 <- ll - f.fff < 1.96
          }
          dx <- -dx
          tol <- F
          while (!tol) {
            dx <- dx / 2
            fff <- fff + dx
            f.fff <- dtlsolver(c(fff, e), t = t, n00 = co[, 3], n10 = co[, 2],
                               n01 = co[, 1], n11 = co[, 4])
            lld <- ll - f.fff
            if (lld < 1.96) dx <- abs(dx) else dx <- -abs(dx)
            if (n > 1000) {
              fff <- NA
              break
            }
            tol <- abs(1.96 - lld) < 1e-3
          }
          cup <- fff

          fff <- upper[2]
          f.fff <- dtlsolver(c(c, fff), t = t, n00 = co[, 3], n10 = co[, 2],
                             n01 = co[, 1], n11 = co[, 4])
          dx <- fff - e
          n <- 0
          lld1.96 <- ll - f.fff < 1.96
          while (lld1.96) {
            fff <- fff + dx
            f.fff <- dtlsolver(c(c, fff), t = t, n00 = co[, 3], n10 = co[, 2],
                               n01 = co[, 1], n11 = co[, 4])
            n <- n + 1
            if (n > 1000) {
              break
            }
            lld1.96 <- ll - f.fff < 1.96
          }
          dx <- -dx
          tol <- F
          while (!tol) {
            dx <- dx / 2
            fff <- fff + dx
            f.fff <- dtlsolver(c(c, fff), t = t, n00 = co[, 3], n10 = co[, 2],
                               n01 = co[, 1], n11 = co[, 4])
            lld <- ll - f.fff
            if (lld < 1.96) dx <- abs(dx) else dx <- -abs(dx)
            if (n > 1000) {
              fff <- NA
              break
            }
            tol <- abs(1.96 - lld) < 1e-3
          }
          eup <- fff
          break
        }
      }
      if (sum(is.na(c(clo, cup, elo, eup))) > 0) {
        warning(paste0("Confidence intervals couldn't be calculated for group ",
                       group,
                       ". Please exclude it from the analysis or use argument 'step' to get an estimate."),
                call. = F)
      }
    } else {

        lli <- ll

        cup <- NULL
        clo <- NULL
        eup <- NULL
        elo <- NULL
        i <- 1
        while ((ll - lli) < 2.1) {
          cup <- c + (step * i)
          lli <- dtlsolver(c(cup, e), t, n00 = co[, 3], n10 = co[, 2],
                           n01 = co[, 1], n11 = co[, 4])
          i   <- i + 1
          if ((ll - lli) < 1.96) next
          break
        }
        i <- 1
        lli <- ll
        while ((ll - lli) < 2.1) {
          clo <- c - (step * i)
          i <- i + 1
          lli <- dtlsolver(c(clo, e), t, n00 = co[, 3], n10 = co[, 2],
                           n01 = co[, 1], n11 = co[, 4])
          if ((ll - lli) < 1.96) next
          break
        }
        i <- 1
        lli <- ll
        while ((ll - lli) < 2.1) {
          eup <- e + (step * i)
          lli <- dtlsolver(c(c, eup), t, n00 = co[, 3], n10 = co[, 2],
                           n01 = co[, 1], n11 = co[, 4])
          i <- i + 1
          if ((ll - lli) < 1.96) next
          break
        }
        i <- 1
        lli <- ll
        while ((ll - lli) < 2.1) {
          elo <- e - (step * i)
          lli <- dtlsolver(c(c, elo), t, n00 = co[, 3], n10 = co[, 2],
                           n01 = co[, 1], n11 = co[, 4])
          i <- i + 1
          if ((ll - lli) < 1.96) next
          break
        }
    }
    sol <- data.frame(c, cup, clo, e, eup, elo, nrow(dataset), -ll)
  } else {
    sol <- data.frame(c, NA, NA, e, NA, NA, nrow(dataset), -ll)
  }
  colnames(sol) <- c("c", "c_up", "c_low", "e", "e_up", "e_low", "N", "NLL")
  sol
}

wrapper2 <- function(list, vectorlist, c, e, step, assembly, jacobian, verbose, CI,
                     group = "NA") {

  ### This function wraps the whole process of obtaining c and e from a
  ### list of datasets with discontinuous times, that could be different for
  ### each dataset.

  time <- vector("list", length(list))
  abs  <- vector("list", length(list))
  pre  <- vector("list", length(list))
  col  <- vector("list", length(list))
  ext  <- vector("list", length(list))
  rows <- 0

  for (i in 1:length(list)) {
    dataset <- list[[i]]
    vector <- vectorlist[[i]]
    co <- data.frame()
    co <- incounts(dataset, vector, assembly)
    time[i] <- data.frame(times(dataset, vector, assembly))
    abs[i]  <- data.frame(co[, 3])
    pre[i]  <- data.frame(co[, 4])
    col[i]  <- data.frame(co[, 2])
    ext[i]  <- data.frame(co[, 1])
    rows <- rows + nrow(dataset)
  }
  if(jacobian){
    if (!requireNamespace("rootSolve", quietly = TRUE)) {
      stop("Package rootSolve needed for this function to work. Please install it.",
           call. = FALSE)
    }
    r <- rootSolve::multiroot(unequal_jacobian_multi, c(c, e), abs = abs,
                              ext = ext, col = col, pre = pre, dt = time,
                              positive = T)
    if (is.nan(r$estim.precis)) {
      stop(paste("Numerical solver didn't converge for", group,
                 " group. Please use lower priors for c and e or exclude it from the analysis."),
           call. = F)
    }
    if (r$estim.precis > 0.0001) {
      warning(paste0("Numerical solver probably converged to a local optimum in group ",
                     group,
                     ". Please use lower priors for c and e or exclude it from the analysis."),
              call. = F)
    }

    c <- r$root[1]
    e <- r$root[2]
    ll <- dtlsolver2(c(c, e), t = time, abs = abs, col = col, ext = ext,
                     pre = pre)
  } else {
    r <- stats::optim(c(c, e), dtlsolver2, t = time, abs = abs, col = col,
                      ext = ext, pre = pre,
                      control = list(maxit = 5000000, fnscale = -1),
                      hessian = F)

    c <- r$par[1]
    e <- r$par[2]
    ll <- r$value
  }
  if (verbose) print(r)
  if (CI) {
    if (is.null(step)) {
      while (CI) {
        mhh <- -matrix(c(hessian_c_m(abs = abs, pre = pre, col = col,
                                     ext = ext, dt = time, c = c, e = e),
                         0,
                         0,
                         hessian_e_m(abs = abs, pre = pre, col = col,
                                        ext = ext, dt = time, c = c, e = e)),
                       nrow = 2)

        h <- suppressWarnings(try(sqrt(diag(solve(mhh))) * 1.96, silent = T))
        if (inherits(h, "try-error")) {
          clo <- NA
          cup <- NA
          elo <- NA
          eup <- NA
          break
        } else {

          if (is.nan(h[1])) h[1] <- .99 * c
          if (is.nan(h[2])) h[2] <- .99 * e

          lower <- c(c - h[1], e - h[2])
          upper <- c(c + h[1], e + h[2])
          n <- 0
          if (lower[1] < 0) {
            fff <- c
            dx <- c
          } else {
            fff <- lower[1]
            f.fff <- dtlsolver2(c(fff, e), t = time, pre=pre, abs = abs,
                                col = col, ext = ext)
            dx <- fff - c
            lld1.96 <- ll - f.fff < 1.96
            while (lld1.96) {
              fff <- fff + dx
              if (fff < 0) {
                fff <- fff - dx
                dx <- fff / 2
                fff <- fff - dx
              }
              f.fff <- dtlsolver2(c(fff, e), t = time, pre=pre, abs = abs,
                                  col = col, ext = ext)
              n <- n + 1
              if (n > 1000) {
                break
              }
              lld1.96 <- ll - f.fff < 1.96
            }
          }
          dx <- -dx
          tol <- F
          while (!tol) {
            dx <- dx / 2
            fff <- fff + dx
            if (fff < 0) {
              fff <- fff - dx
              dx <- fff / 2
              fff <- fff - dx
            }
            f.fff <- dtlsolver2(c(fff, e), t = time, pre = pre, abs = abs,
                                col = col, ext = ext)
            lld <- ll - f.fff
            if (lld < 1.96) dx <- -abs(dx) else dx <- abs(dx)
            n <- n + 1
            if (n > 1000) {
              fff <- NA
              break
            }
            tol <- abs(1.96 - lld) < 1e-3
          }
          clo <- fff

          n <- 0
          if (lower[2] < 0) {
            fff <- e
            dx <- e
          } else {
            fff <- lower[2]
            f.fff <- dtlsolver2(c(c, fff), t = time, pre = pre, abs = abs,
                                col = col, ext = ext)
            dx <- fff - e
            lld1.96 <- ll - f.fff < 1.96
            while (lld1.96) {
              fff <- fff + dx
              if (fff < 0) {
                fff <- fff - dx
                dx <- fff / 2
                fff <- fff - dx
              }
              f.fff <- dtlsolver2(c(c, fff), t = time, pre = pre, abs = abs,
                                  col = col, ext = ext)
              n <- n + 1
              if (n > 1000) {
                break
              }
              lld1.96 <- ll - f.fff < 1.96
            }
          }
          dx <- -dx
          tol <- F
          while (!tol) {
            dx <- dx / 2
            fff <- fff + dx
            if (fff < 0) {
              fff <- fff - dx
              dx <- fff / 2
              fff <- fff - dx
            }
            f.fff <- dtlsolver2(c(c, fff), t = time, pre = pre, abs = abs,
                                col = col, ext = ext)
            lld <- ll - f.fff
            if (lld < 1.96) dx <- -abs(dx) else dx <- abs(dx)
            if (n > 1000) {
              fff <- NA
              break
            }
            tol <- abs(1.96 - lld) < 1e-3
          }
          elo <- fff

          n <- 0
          fff <- upper[1]
          f.fff <- dtlsolver2(c(fff, e), t = time, pre = pre, abs = abs,
                              col = col, ext = ext)
          dx <- fff - c
          lld1.96 <- ll - f.fff < 1.96
          while (lld1.96) {
            fff <- fff + dx
            f.fff <- dtlsolver2(c(fff, e), t = time, pre = pre, abs = abs,
                                col = col, ext = ext)
            n <- n + 1
            if (n > 1000) {
              break
            }
            lld1.96 <- ll - f.fff < 1.96
          }
          dx <- -dx
          tol <- F
          while (!tol) {
            dx <- dx / 2
            fff <- fff + dx
            f.fff <- dtlsolver2(c(fff, e), t = time, pre = pre, abs = abs,
                                col = col, ext = ext)
            lld <- ll - f.fff
            if (lld < 1.96) dx <- abs(dx) else dx <- -abs(dx)
            if (n > 1000) {
              fff <- NA
              break
            }
            tol <- abs(1.96 - lld) < 1e-3
          }
          cup <- fff

          fff <- upper[2]
          f.fff <- dtlsolver2(c(c, fff), t = time, pre = pre, abs = abs,
                              col = col, ext = ext)
          dx <- fff - e
          n <- 0
          lld1.96 <- ll - f.fff < 1.96
          while (lld1.96) {
            fff <- fff + dx
            f.fff <- dtlsolver2(c(c, fff), t = time, pre = pre, abs = abs,
                                col = col, ext = ext)
            n <- n + 1
            if (n > 1000) {
              break
            }
            lld1.96 <- ll - f.fff < 1.96
          }
          dx <- -dx
          tol <- F
          while (!tol) {
            dx <- dx / 2
            fff <- fff + dx
            f.fff <- dtlsolver2(c(c, fff), t = time, pre = pre, abs = abs,
                                col = col, ext = ext)
            lld <- ll - f.fff
            if (lld < 1.96) dx <- abs(dx) else dx <- -abs(dx)
            if (n > 1000) {
              fff <- NA
              break
            }
            tol <- abs(1.96 - lld) < 1e-3
          }
          eup <- fff
          break
        }
      }
      if (sum(is.na(c(clo, cup, elo, eup))) > 0) {
        warning(paste0("Confidence intervals couldn't be calculated for group ",
                       group,
                       ". Please exclude it from the analysis or use argument 'step' to get an estimate."),
                call. = F)
      }
    } else {

      lli <- ll

      cup <- NULL
      clo <- NULL
      eup <- NULL
      elo <- NULL
      i <- 1
      while ((ll - lli) < 2.1) {
        cup <- c + (step * i)
        lli <- dtlsolver2(c(cup, e), t = time, pre = pre, abs = abs, col = col,
                          ext = ext)
        i   <- i + 1
        if ((ll - lli) < 1.96) next
        break
      }
      i <- 1
      lli <- ll
      while ((ll - lli) < 2.1) {
        clo <- c - (step * i)
        i <- i + 1
        lli <- dtlsolver2(c(clo, e), t = time, pre = pre, abs = abs, col = col,
                          ext = ext)
        if ((ll - lli) < 1.96) next
        break
      }
      i <- 1
      lli <- ll
      while ((ll - lli) < 2.1) {
        eup <- e + (step * i)
        lli <- dtlsolver2(c(c, eup), t = time, pre = pre, abs = abs, col = col,
                          ext = ext)
        i <- i + 1
        if ((ll - lli) < 1.96) next

        break
      }
      i <- 1
      lli <- ll
      while ((ll - lli) < 2.1) {
        elo <- e - (step * i)
        lli <- dtlsolver2(c(c, elo), t = time, pre = pre, abs = abs, col = col,
                          ext = ext)
        i <- i + 1
        if ((ll - lli) < 1.96) next
        break
      }
    }
    sol <- data.frame(c, cup, clo, e, eup, elo, rows, -ll)
  } else {
    sol <- data.frame(c, NA, NA, e, NA, NA, rows, -ll)
  }
  colnames(sol) <- c("c", "c_up", "c_low", "e", "e_up", "e_low", "N", "NLL")
  sol
}

unequal_jacobian <- function(x, N00, N01, N10, N11, dt){
  c <- x[1]
  e <- x[2]
  F1 <- (e * (N10 + N11) - c * (N00 + N01))/(c * (c + e)) +
    (dt * (N01 + N10)) / (-1 + exp(dt * (c + e))) +
    (N00 - c * N00 * dt)/(c + e * exp(dt * (c + e))) -
    (e * N11 * (1 + c * dt))/(c * (e + c * exp(dt * (c + e))))

  F2 <- (-e * (N10 + N11) + c * (N00 + N01))/(e * (c + e)) +
    (dt * (N01 + N10)) / (-1 + exp(dt * (c + e))) +
    (N11 - e * N11 * dt)/(e + c * exp(dt * (c + e))) -
    (c * N00 * (1 + e * dt))/(e * (c + e * exp(dt * (c + e))))
  c(F1 = sum(F1), F2 = sum(F2))
}


unequal_jacobian_multi <- function(x, abs, ext, col, pre, dt){
  c <- x[1]
  e <- x[2]
  F1 <- 0
  F2 <- 0
  for(i in 1:length(dt)){
  temp1 <- (e * (col[[i]] + pre[[i]]) - c * (abs[[i]] + ext[[i]]))/(c *
                                                                      (c + e)) +
    (dt[[i]] * (ext[[i]] + col[[i]])) / (-1 + exp(dt[[i]] * (c + e))) +
    (abs[[i]] - c * abs[[i]] * dt[[i]])/(c + e * exp(dt[[i]] * (c + e))) -
    (e * pre[[i]] * (1 + c * dt[[i]]))/(c * (e + c * exp(dt[[i]] * (c + e))))

  temp2 <- (-e * (col[[i]] + pre[[i]]) + c * (abs[[i]] + ext[[i]]))/(e *
                                                                      (c + e)) +
    (dt[[i]] * (ext[[i]] + col[[i]])) / (-1 + exp(dt[[i]] * (c + e))) +
    (pre[[i]] - e * pre[[i]] * dt[[i]])/(e + c * exp(dt[[i]] * (c + e))) -
    (c * abs[[i]] * (1 + e * dt[[i]]))/(e * (c + e * exp(dt[[i]] * (c + e))))
  F1 <- F1 + sum(temp1)
  F2 <- F2 + sum(temp2)
  }
  c(F1 = F1, F2 = F2)
}


hessian_c_m <- function(abs, pre, ext, col, dt, c, e){
  hc <- 0

  for(i in 1:length(dt)){
  Edt <- exp(dt[[i]] * (c + e))

  temphc <- ext[[i]] / (c + e)^2 -
    (col[[i]] * e) / (c * (c + e)^2) -
    (pre[[i]] * e) / (c * (c + e)^2) -
    (col[[i]] * e) / (c^2 * (c + e)) -
    (pre[[i]] * e) / (c^2 * (c + e)) +
    (pre[[i]] * e) / (c^2 * (Edt * c + e)) +
    abs[[i]] * (1 / (c + e)^2 - 1 / (c + Edt * e)^2 -
                  (Edt * dt[[i]] * e) / (c + Edt * e)^2) +
    Edt * (-((ext[[i]] * dt[[i]]^2) / (-1 + Edt)^2) -
             (col[[i]] * dt[[i]]^2) / (-1 + Edt)^2 +
             (pre[[i]] * (1 + dt[[i]] * c)^2 * e) / (c * (Edt * c + e)^2) +
             (abs[[i]] * dt[[i]] * (-1 + dt[[i]] * c) * e) / (c + Edt * e)^2)
  hc <- hc + sum(temphc)
  }
  hc
}


hessian_e_m <- function(abs, pre, ext, col, dt, c, e){
  he <- 0
  for(i in 1:length(dt)){
  Edt <- exp(dt[[i]] * (c + e))

  temphe <- col[[i]] / (c + e)^2 +
    pre[[i]] / (c + e)^2 -
    (abs[[i]] * c) / (e * (c + e)^2) -
    (ext[[i]] * c) / (e * (c + e)^2) -
    (abs[[i]] * c) / (e^2 * (c + e)) -
    (ext[[i]] * c) / (e^2 * (c + e)) -
    pre[[i]] / (Edt * c + e)^2 +
    (pre[[i]] * dt[[i]] * e) / (Edt * c + e)^2 -
    (pre[[i]] * dt[[i]]) / (Edt * c + e) +
    (abs[[i]] * c) / (e^2 * (c + Edt * e)) +
    Edt * (-((ext[[i]] * dt[[i]]^2) / (-1 + Edt)^2) -
             (col[[i]] * dt[[i]]^2) / (-1 + Edt)^2 +
             (pre[[i]] * dt[[i]] * c * (-1 + dt[[i]] * e)) / (Edt * c + e)^2 +
             (abs[[i]] * c * (1 + dt[[i]] * e)^2) / (e * (c + Edt * e)^2))
  he <- he + sum(temphe)
  }
  he
}
