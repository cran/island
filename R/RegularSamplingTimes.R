#' c/e rates for a regular sampling scheme
#'
#' \code{regular_sampling_scheme} estimates colonization and extinction rates
#' for a community or groups in a community. \cr \code{NLL_rss} returns the Negative
#' Log-Likelihood of a pair of colonization and extinction rates for a regular
#' sampling scheme.
#'
#' @param x A single dataset.
#' @param vector A vector indicating the columns with presence-absence data.
#' @param level The name of the column that contain groups used to subset them
#'   and calculate their colonization and extinction rates.
#' @param n Minimal number of rows for each group.
#' @param step Accuracy to calculate the c_e pairs with.
#' @param CI Logical. Should confidence intervals be returned?
#' @param c A colonization rate.
#' @param e An extinction rate.
#' @examples regular_sampling_scheme(alonso15[[1]], 3:6)
#' regular_sampling_scheme(alonso15[[1]], 3:6, "Guild", n = 5)
#' regular_sampling_scheme(alonso15[[1]], 3:6, "Guild", n = 5, CI = TRUE)
#' NLL_rss(alonso15[[1]], 3:6, 0.52, 0.39)
#' @return \code{regular_sampling_scheme} returns a dataframe with colonization and extinction rates along with their
#'   lower and upper confidence intervals (optional), for each group if
#'   specified, and its number of rows and NLL.
#'   \code{NLL_rss} gives the NLL for a dataframe given a specific c and e.
#' @details The confidence intervals are calculated with a binary search seeded
#'   with the hessian of the estimated rates.
#' @seealso \code{\link{irregular_single_dataset}},
#'   \code{\link{irregular_multiple_datasets}}
#' @export
regular_sampling_scheme <- function(x, vector, level = NULL, n = NULL,
                                    step = NULL, CI = F) {

  if (is.null(level)) {
    rates(changes(x,vector), CI)
  } else {
    if (is.null(step)) {
      rlevel2(x, level, vector, n, CI)
    } else {
      rlinterval(x, level, vector, n, step, CI)
    }
  }
}

#' @rdname regular_sampling_scheme
#' @export
NLL_rss <- function(x, vector, c, e) {
  -lfun2(c, e, changes(x, vector))
}

changes <- function(x, vector) {

  ### This function identifies and counts the kind of transition produced in a
  ### dataframe "x" among every column in "vector". The output is a dataframe
  ### with the number of extinctions, colonizations and mantained absences or
  ### presences in the dataframe.

  N01 <- 0
  N10 <- 0
  N00 <- 0
  N11 <- 0
  resultado <- data.frame()
  for (i in 1:nrow(x)){
    for (j in 1:(length(vector)-1)){
      if (utils::tail(vector, 1) == j) break
      if (x[i, vector[j]] < x[i, vector[j + 1]]) N10 <- N10 + 1
      if (x[i, vector[j]] > x[i, vector[j + 1]]) N01 <- N01 + 1
      if (x[i, vector[j]] == x[i, vector[j + 1]] && x[i, vector[j]] == 0) N00 <- N00 + 1
      if (x[i, vector[j]] == x[i, vector[j + 1]] && x[i, vector[j]] == 1) N11 <- N11 + 1
      if (i != nrow(x)) next
      resultado <- rbind(c(N01, N10, N00, N11, nrow(x)))
    }
  }
  resultado
}

lfun <- function(c, e, x) {

  ### This function calculates the likelihood of a given set of transitions
  ### x for the rates c and e.

  (1 - ( (c / (e+c)) * (1 - exp( - (e + c))))) ^ x[1, 3] * ( (c / (e + c)) *
  (1 - exp( - (e + c)))) ^ x[1, 2] * ( (e / (e + c)) * (1 - exp( - (e + c)))) ^
  x[1, 1] * (1 - ( (e / (e + c)) * (1 - exp( - (e + c))))) ^ x[1, 4]
}

lfun2 <- function(c, e, x) {

  ### This function calculates the log-likelihood of a given set of
  ### transitions x for the rates c and e.

  log(1 - ( (c / (e + c)) * (1 - exp( - (e + c))))) * x[1, 3] +
  log(1 - ( (e / (e + c)) * (1 - exp( - (e + c))))) * x[1, 4] +
  log( (c / (e + c)) * (1 - exp( - (e + c)))) * x[1, 2] +
  log( (e / (e + c)) * (1 - exp( - (e + c)))) * x[1, 1]

}

lfunr <- function(T01, T10, x) {

  ### This function calculates the log-likelihood of a given set of
  ### transitions x for the transition probabilities T01 and T10.

  log(1 - T10) * x[1, 3] + log(T10) * x[1, 2] + log(T01) * x[1, 1] +
  log(1 - T01) * x[1, 4]
}

lfunrsolver <- function(y, x) {

  ### This function calculates the loglikelihood of a given dataset of
  ### transitions x for the vector y of transition probabilities T01
  ### and T10.

  T01 <- y[1]
  T10 <- y[2]
  log(1 - T10) * x[1, 3] + log(T10) * x[1, 2] + log(T01) * x[1, 1] +
  log(1 - T01) * x[1, 4]
}

like <- function(x) {

  ### This function calculates the log-likelihood of a given set of
  ### transitions x for the c and e calculated with the transitions in
  ### x.

  out <- numeric()
  edivc1  <- x[1, 1] * (x[1, 2] + x[1, 3]) / (x[1, 2] * (x[1, 1] + x[1, 4]))
  eplusc1 <- - log( (x[1, 3] * x[1, 4] - x[1, 2] * x[1, 1]) /
                    ( (x[1, 2] + x[1, 3]) * (x[1, 1] + x[1, 4])))
  c <- eplusc1 / (edivc1 + 1)
  e <- eplusc1 / (1 + (1 / edivc1))
  T00 <- 1 - ( (c / (e + c)) * (1 - exp( - (e + c))))
  T11 <- 1 - ( (e / (e + c)) * (1 - exp( - (e + c))))
  T01 <- (e / (e + c)) * (1 - exp( - (e + c)))
  T10 <- (c / (e + c)) * (1 - exp( - (e + c)))
  out <- log(T00) * x[1, 3] + log(T10) * x[1, 2] + log(T01) * x[1, 1] +
         log(T11) * x[1, 4]
  out
}

rates <- function(y, CI) {

  ### This function calculates the rates and probabilities (c, e, T00, T01,
  ### T11, T10) for the given set of transitions x.

  sol <- data.frame()
  edivc1 <- (y[1, 1] * (y[1, 2] + y[1, 3])) / (y[1, 2] * (y[1, 1] + y[1, 4]))
  eplusc1 <- - log( ( (y[1, 3] * y[1, 4]) - (y[1, 2] * y[1, 1])) /
                    ( (y[1, 2] + y[1, 3]) * (y[1, 1] + y[1, 4])))
  c <- eplusc1 / (edivc1 + 1)
  e <- eplusc1 / (1 + (1 / edivc1))

  ll <- lfun2(c, e, x = y)
  if(CI){
    h <- sqrt(diag(solve(-matrix(c(hessian_c_r(n00 = y[3], n11 = y[4],
                              n01 = y[1], n10 = y[2], dt = 1, c = c, e = e), 0,
                              0, hessian_e_r(n00 = y[3], n11 = y[4], n01 = y[1],
                              n10 = y[2], dt = 1, c = c, e = e)), nrow = 2)))) *
      1.96
    lower <- c(c - h[1], e - h[2])
    upper <- c(c + h[1], e + h[2])

    if(lower[1]<0){
      fff <- lower[1]
      dx <- fff - c

    } else {
      fff <- lower[1]
      f.fff <- lfun2(fff, e, x = y)

      dx <- fff - c

      lld1.96 <- ll - f.fff < 1.96

      while(lld1.96){

        fff <- fff + dx
        f.fff <- lfun2(fff, e, x = y)
        lld1.96 <- ll - f.fff < 1.96
      }
    }
    dx <- -dx
    tol <- F
    while(!tol){
      dx <- dx / 2
      fff <- fff + dx
      if(fff < 0){
        fff <- fff - dx
        dx <- fff / 2
        fff <- fff - dx
      }
      f.fff <- lfun2(fff, e, x = y)
      lld <- ll - f.fff
      if (lld < 1.96) dx <- -abs(dx) else dx <- abs(dx)
      tol <- abs(1.96 - lld) < 1e-3
    }
    clo <- fff

    if(lower[2]<0){
      fff <- lower[2]
      dx <- fff - e

    } else {
      fff <- lower[2]
      f.fff <- lfun2(c, fff, x = y)

      dx <- fff - e

      lld1.96 <- ll - f.fff < 1.96

      while(lld1.96){

        fff <- fff + dx
        f.fff <- lfun2(c, fff, x = y)
        lld1.96 <- ll - f.fff < 1.96
      }
    }
    dx <- -dx
    tol <- F
    while(!tol){
      dx <- dx / 2
      fff <- fff + dx
      if(fff < 0){
        fff <- fff - dx
        dx <- fff / 2
        fff <- fff - dx
      }
      f.fff <- lfun2(c, fff, x = y)
      lld <- ll - f.fff
      if (lld < 1.96) dx <- -abs(dx) else dx <- abs(dx)
      tol <- abs(1.96 - lld) < 1e-3
    }
    elo <- fff

    fff <- upper[1]

    f.fff <- lfun2(fff, e, x = y)
    dx <- fff - c

    lld1.96 <- ll - f.fff < 1.96
    while(lld1.96){

      fff <- fff + dx
      f.fff <- lfun2(fff, e, x = y)
      lld1.96 <- ll - f.fff < 1.96
    }
    dx <- -dx
    tol <- F
    while(!tol){
      dx <- dx / 2
      fff <- fff + dx
      f.fff <- lfun2(fff, e, x = y)
      lld <- ll - f.fff
      if (lld < 1.96) dx <- abs(dx) else dx <- -abs(dx)
      tol <- abs(1.96 - lld) < 1e-3
    }
    cup <- fff


    fff <- upper[2]

    f.fff <- lfun2(c, fff, x = y)
    dx <- fff - e

    lld1.96 <- ll - f.fff < 1.96
    while(lld1.96){

      fff <- fff + dx
      f.fff <- lfun2(c, fff, x = y)
      lld1.96 <- ll - f.fff < 1.96
    }
    dx <- -dx
    tol <- F
    while(!tol){
      dx <- dx / 2
      fff <- fff + dx
      f.fff <- lfun2(c, fff, x = y)
      lld <- ll - f.fff
      if (lld < 1.96) dx <- abs(dx) else dx <- -abs(dx)
      tol <- abs(1.96 - lld) < 1e-3
    }
    eup <- fff
    sol <- data.frame(c, cup, clo, e, eup, elo, y[5], -ll)
  } else {
    sol <- data.frame(c, NA, NA, e, NA, NA, y[5], -ll)
  }

  colnames(sol) <- c("c", "c_up", "c_low", "e", "e_up", "e_low", "N", "NLL")
  sol

}

rlevel2 <- function(x, level, vector, n, CI) {

  ### Starting from a dataset x that contains transitions in the columns
  ### specified in vector, this function calculates the rates and
  ### probabilities (c, e, T00, T01, T11, T10)  for the groups in level
  ### that have more than a number otus of replicas.

  results <- data.frame()
  ## First, we make an index of the actual elements in "level".
  index <- levels(factor(x[, level]))
  ## Now, start the calculations.
  for (i in 1:length(index)) {
    y <- x[x[, level] == index[i], ]
    if (nrow(y) < n) next
    rr <- rates(changes(y, vector), CI)
    rr <- data.frame(index[i], rr)
    results <- rbind(results, rr)
  }
  colnames(results) <- c("Group", "c", "c_up", "c_low", "e", "e_up", "e_low", "N", "NLL")
  results
}

rlinterval <- function(x, level, vector, n, step, CI) {

  ### Starting from a dataset x that contains transitions in the columns
  ### specified in vector, this function calculates the rates and its 95%
  ### confidence interval for the groups in level that have more than a
  ### number "otus" of replicas.

  out <- data.frame()
  out2 <- data.frame()
  ## First, we make an index of the actual elements in "level".
  index <- levels(factor(x[, level]))
  ## Now, start the calculations.
  for (i in 1:length(index)) {
    y <- x[x[, level] == index[i], ]
    if (nrow(y) < n) next
    rr <- rates(changes(y, vector), CI)
    c1 <- rr[1]
    e1 <- rr[4]
    if (is.na(c1) == T) next
    z <- changes(y, vector)

    llce <- lfun2(c1, e1, z)
    df <- data.frame(matrix(0, 1, 6))
    lli <- llce

    cup <- NULL
    clo <- NULL
    eup <- NULL
    elo <- NULL

    j <- 1
    while ( (llce-lli) < 2) {
      cup <- c1 + (step * j)
      lli <- lfun2(cup, e1, z)
      j <- j + 1
      if ( (llce - lli) < 1.96) next
      df[1, 3] <- cup
      break
    }
    j <- 1
    lli <- llce
    while ( (llce - lli) < 2) {
      clo <- c1 - (step * j)
      j <- j + 1
      lli <- lfun2(clo, e1, z)
      if ( (llce - lli) < 1.96) next
      df[1, 2] <- clo
      break
    }
    j <- 1
    lli <- llce
    while ( (llce - lli) < 2) {
      eup <- e1 + (step * j)
      lli <- lfun2(c1, eup, z)
      j <- j + 1
      if ( (llce - lli) < 1.96) next
      df[1, 6] <- eup
      break
    }
    j <- 1
    lli <- llce
    while ( (llce - lli) < 2) {
      elo <- e1 - (step * j)
      lli <- lfun2(c1, elo, z)
      j <- j + 1
      if ( (llce - lli) < 1.96) next
      df[1, 5] <- elo
      break
    }

    out <- data.frame(index[i], c1, df[1, 3], df[1, 2], e1, df[1, 6],
                      df[1, 5], nrow(y), rr[8])
    colnames(out) <- c("Group", "c", "c_up", "c_low", "e", "e_up", "e_low", "N", "NLL")
    out2 <- rbind(out2, out)


  }
  out2
}


hessian_c_r <- function(n00, n11, n01, n10, dt, c, e){
  Edt <- exp(dt * (c + e))

  n01/(c + e)^2 -
  (n10 * e)/(c * (c + e)^2) -
  (n11 * e)/(c * (c + e)^2) -
  (n10 * e)/(c^2 * (c + e)) -
  (n11 * e)/(c^2 * (c + e)) +
  (n11 * e)/(c^2 * (Edt * c + e)) +
  n00 * (1/(c + e)^2 - 1/(c + Edt * e)^2 - (Edt * dt * e)/(c + Edt * e)^2) +
  Edt * (-((n01 * dt^2)/(-1 + Edt)^2) - (n10 * dt^2)/(-1 + Edt)^2 + (n11 * (1 +
    dt * c)^2 * e)/(c * (Edt * c + e)^2) + (n00 * dt * (-1 + dt * c) * e)/(c +
    Edt * e)^2)
}

hessian_e_r <- function(n00, n11, n01, n10, dt, c, e){
  Edt <- exp(dt * (c + e))

  n10/(c + e)^2 +
  n11/(c + e)^2 -
  (n00 * c)/(e * (c + e)^2) -
  (n01 * c)/(e * (c + e)^2) -
  (n00 * c)/(e^2 * (c + e)) -
  (n01 * c)/(e^2 * (c + e)) -
  n11/(Edt * c + e)^2 +
  (n11 * dt * e)/(Edt * c + e)^2 -
  (n11 * dt)/(Edt * c + e) +
  (n00 * c)/(e^2 * (c + Edt * e)) +
  Edt * (-((n01 * dt^2)/(-1 + Edt)^2) - (n10 * dt^2)/(-1 + Edt)^2 + (n11 * dt *
    c * (-1 + dt * e))/(Edt * c + e)^2 + (n00 * c * (1 + dt * e)^2)/(e * (c +
    Edt * e)^2))

}
