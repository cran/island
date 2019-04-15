context("Detectability")

### The next test works locally but is not working when checking.

# library(unmarked)
#
# test_that("island produce the same estimates as unmarked when presented with the same data",
#           {
#             Data1 <- lakshadweep[[1]][, 1:16]
#             Name_of_Factors <- c("Species","Atoll","Guild")
#             Factors <- Filter(is.factor, Data1)
#             No_of_Factors <- length(Factors[1,])
#             n <- No_of_Factors + 1
#             D1 <- as.matrix(Data1[1:nrow(Data1),n:ncol(Data1)])
#             Time <- as.double(D1[1,])
#             P1 <- as.matrix(D1[2:nrow(D1),1:ncol(D1)])
#             # Dealing with time.
#             Time_Vector <- as.numeric(names(table(Time)))
#             Transects   <- as.numeric((table(Time)))
#             R1 <- sss_cedp(P1, Time_Vector, Transects,
#                            Colonization=0.5, Extinction=0.5, Detectability=0.5,
#                            Phi_Time_0=0.5,
#                            Tol=1.0e-8, Verbose = F)
#
#             # Preparing data for function unmarked::colext
#             colNA <- rep(NA, 156)
#             yy <- cbind(P1[, 1:2], colNA, colNA, P1[, 3:5], colNA, P1[, 6:13]) #
#             year <- matrix(Time_Vector,
#                            nrow(yy), 4, byrow=TRUE) #Indicates which years the sites have been surveyed
#
#             simUMF <- unmarkedMultFrame(
#               y = yy,
#               yearlySiteCovs = list(year = year),
#               numPrimary=4)
#
#             m0 <- colext(psiformula= ~1, gammaformula = ~ 1, epsilonformula = ~ 1,
#                          pformula = ~ 1, data = simUMF, method="BFGS")
#             summary(m0)
#             expect_lt(R1$P - plogis(m0@estimates@estimates$psi@estimates),
#                       .05 * R1$P)
#             expect_lt(R1$D - plogis(m0@estimates@estimates$det@estimates),
#                       .05 * R1$D)
#             tp <-  cetotrans(R1$C, R1$E)
#             expect_lt(tp[1] - plogis(m0@estimates@estimates$ext@estimates),
#                       .05 * tp[1])
#             expect_lt(tp[2] - plogis(m0@estimates@estimates$col@estimates),
#                       .05 * tp[2])
#             expect_lt(abs(R1$NLL - m0@negLogLike), 1)
#           })
#
# detach("package:unmarked")

context("Testing detectability")

test_that("the detectability functions produce the same estimates", {
  Data <- lakshadweepPLUS[[1]]
  Guild_Tag = c("Alg","Cor","Mac","Mic","Omn","Pis","Zoo")
  Time <- (as.vector(c(2000, 2000, 2001, 2001, 2001, 2001, 2002, 2002, 2002, 2002, 2003, 2003, 2003, 2003, 2010, 2010, 2011, 2011, 2011, 2011, 2012, 2012, 2013, 2013, 2013, 2013)))
  Data <- Data[Data$Atoll == "KADMATH", -c(26, 27)]
  Rm <- mss_cedp(Data, Time, Factor = 3, Tags = Guild_Tag,
                 PerfectDetectability = FALSE, z = 4, Verbose = 0)
  Ru <- upgma_model_selection(Data, Time, Factor = 3, Tags = Guild_Tag,
                              PerfectDetectability = FALSE, z = 4, Verbose = 0)

  NLL <- sum(Rm[, 5]) # This should be equivalent to the NLL of the upgma's last model

  expect_equal(NLL, Ru[7, 2], tolerance = 1)

  DataB <- Data[Data$Guild == "Macroinvertivore", 4:29]
  Time2 <- unique(Time)
  Transects   <- c(2, 4, 4, 4, 2, 4, 2, 4)
  Rs <- sss_cedp(DataB, Time2, Transects,
                 Colonization=0.5, Extinction=0.5, Detectability=0.5,
                 Phi_Time_0=0.5,
                 Tol=1.0e-12, Verbose = 0)

  expect_equivalent(Rs$NLL, Rm[3, 5]) # The NLL for the macroinvertivores is the same with mss and sss.
})
context("Name of test context")
