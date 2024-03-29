% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/da_MacKenzie_mle.R
\name{sss_cedp}
\alias{sss_cedp}
\title{MacKenzie etal (2003) likelihood approach for estimating colonization/extinction
parameters (with imperfect detectability)}
\usage{
sss_cedp(
  Data,
  Time,
  Transects,
  Colonization = 0.1,
  Extinction = 0.1,
  Detectability = 0.99,
  Phi_Time_0 = 0.5,
  Tol = 1e-06,
  MIT = 100,
  C_MAX = 2,
  C_min = 0,
  E_MAX = 2,
  E_min = 0,
  D_MAX = 0.999,
  D_min = 0.001,
  P_MAX = 0.999,
  P_min = 0.001,
  I_0 = 0,
  I_1 = 1,
  I_2 = 2,
  I_3 = 3,
  z = 4,
  Verbose = 0,
  Minimization = TRUE
)
}
\arguments{
\item{Data}{S x N matrix containing presence data per transect (in cols):}

\item{Time}{an array of length n containing increasing sampling times (without repetitions)}

\item{Transects}{an integer array of length n containing the number of transects per sampling time}

\item{Colonization}{guess value to initiate search / parameter value}

\item{Extinction}{guess value to initiate search / parameter value}

\item{Detectability}{guess value to initiate search / param
eter value}

\item{Phi_Time_0}{guess value to initiate search / parameter value}

\item{Tol}{Stopping criteria of the search algorithm}

\item{MIT}{max number of iterations of the search algorithm}

\item{C_MAX}{max value of colonization values}

\item{C_min}{min value of colonization values}

\item{E_MAX}{max value of extinction values}

\item{E_min}{min value of extinction values}

\item{D_MAX}{max value of detectability values}

\item{D_min}{min value of detectability values}

\item{P_MAX}{max value for the initial presence probability on the site}

\item{P_min}{min value for the initial presence probability on the site}

\item{I_0}{parameter index of 1st parameter}

\item{I_1}{parameter index of 2nd parameter}

\item{I_2}{parameter index of 3rd parameter}

\item{I_3}{parameter index of 4th parameter}

\item{z}{dimension of the parameter subspace}

\item{Verbose}{more/less (1/0) output information}

\item{Minimization}{TRUE/FALSE.}
}
\value{
A list with five components (Colonization, Extinction, Detectability,
P_0, and Negative Log-Likelihood).
}
\description{
\code{sss_cedp} conducts a maximum likelihood estimation of model parameters
(Colonization, Extinction, Detectability, and Phi_Time_0) of MacKenzie et al
(2003) colonization-extinction model. This function is an alternative to
\code{mss_cedp} that takes a different input (a 2D array), and requires the same
sampling structure for all input data matrix rows, this is, no missing data
defining a heterogeneous sampling structure across rows are allowed. As an advantage,
it may run faster than \code{mss_cedp}.
}
\details{
Maximum likelihood parameter estimation is conducted through bounded searches.
This is the reason why the minimum and maximum values for each axis should be given
as input arguments. The optimization procedure is the simplex method. A bounded
parameter space implies that in case a neg loglikelihood (NLL) evaluation is
required outside from these boundaries, the returned value for this NLL evaluation
is artifically given as the maximum number the machine can hold.
The array Parameters (I_0, I_1, I_2, I_3) has to be a permutation of (0, 1, 3, 4).
This parameter indeces along with the imput parameter 'z' are used to define a
subparameter space where the search will be conducted. If z = 2, then the search
will take place on the plane defined by model parameters (I_0, I_1). These indeces
are model parameter keys: colonization (0), extinction (1), detectability (2), and
Phi_Time_0 (3). For instance, if (I_0, I_1, I_2, I_3) is (2, 3, 1, 0), and z = 2,
then the search will take place whithin the subparemeter space defined by the
detection probability (Detectability) and the probability of presence at time 0
(Phi_Time_0). If Minimization is TRUE (default value), then the whole mle is
conducted. If FALSE, the function only return the NLL value at the input model
parameter values. Likelihood evaluations are exact provided the number of 'absences'
corresponding to either true absences or undetected presences in the input data
matrix is not to high.
}
\examples{
Data1 <- lakshadweep[[1]]
Name_of_Factors <- c("Species","Atoll","Guild")
Factors <- Filter(is.factor, Data1)
No_of_Factors <- length(Factors[1,])
n <- No_of_Factors + 1
D1 <- as.matrix(Data1[1:nrow(Data1),n:ncol(Data1)])
Time <- as.double(D1[1,])
P1 <- as.matrix(D1[2:nrow(D1),1:ncol(D1)])
# Dealing with time.
Time_Vector <- as.numeric(names(table(Time)))
Transects   <- as.numeric((table(Time)))
R1 <- sss_cedp(P1, Time_Vector, Transects,
                       Colonization=0.5, Extinction=0.5, Detectability=0.5,
                       Phi_Time_0=0.5,
                       Tol=1.0e-8, Verbose = 1)

}
