#' MacKenzie etal (2003) likelihood approach for estimating colonization/extinction
#' parameters (with impecfect detectability) 
#' 
#' \code{da_MacKenzie_mle} conducts a maximum likelihood estimation of model parameters
#' (Colonization, Extinction, Detectability, and Phi_Time_0) of MacKenzie et al
#' (2003) colonization-extinction model. This function is an alternative to
#' \code{da_mle_cedp} that takes a different input (a 2D array), and requires the same
#' sampling structure for all input data matrix rows, this is, no missing data
#' defining a heterogeneous sampling structure across rows are allowed. As an advantage,
#' it may run faster than \code{da_mle_cedp}. 
#'
## Details:
#' Maximum likelihood parameter estimation is conducted through bounded searches.
#' This is the reason why the minimum and maximum values for each axis should be given
#' as input arguments. The optimization procedure is the simplex method. A bounded
#' parameter space implies that in case a neg loglikelihood (NLL) evaluation is
#' required outside from these boundaries, the returned value for this NLL evaluation
#' is artifically given as the maximum number the machine can hold. 
#' The array Parameters (I_0, I_1, I_2, I_3) has to be a permutation of (0, 1, 3, 4).
#' This parameter indeces along with the imput parameter 'z' are used to define a 
#' subparameter space where the search will be conducted. If z = 2, then the search 
#' will take place on the plane defined by model parameters (I_0, I_1). These indeces 
#' are model parameter keys: colonization (0), extinction (1), detectability (2), and
#' Phi_Time_0 (3). For instance, if (I_0, I_1, I_2, I_3) is (2, 3, 1, 0), and z = 2,
#' then the search will take place whithin the subparemeter space defined by the
#' detection probability (Detectability) and the probability of presence at time 0
#' (Phi_Time_0). If Minimization is TRUE (default value), then the whole mle is
#' conducted. If FALSE, the function only return the NLL value at the input model
#' parameter values. Likelihood evaluations are exact provided the number of 'absences'
#' corresponding to either true absences or undetected presences in the input data
#' matrix is not to high.  
#' 
##  Input:
#' @param  Data S x N matrix containing presence data per transect (in cols):
#' @param  Time an array of length n containing sampling times (without repetitions)
#' @param  Transects an integer array of length n containing the number of transects per sampling time 
#' @param  Colonization guess value to initiate search / parameter value
#' @param  Extinction guess value to initiate search / parameter value
#' @param  Detectability guess value to initiate search / param
#' eter value
#' @param  Phi_Time_0 guess value to initiate search / parameter value 
#' @param  Tol Stopping criteria of the search algorithm 
#' @param  MIT max number of iterations of the search algorithm
#' @param  C_min min value of colonization values
#' @param  C_MAX max value of colonization values
#' @param  E_min min value of extinction values
#' @param  E_MAX max value of extinction values
#' @param  D_min min value of detectability values	
#' @param  D_MAX max value of detectability values
#' @param  P_min min value for the initial presence probability on the site
#' @param  P_MAX max value for the initial presence probability on the site  
#' @param  I_0 parameter index of 1st parameter
#' @param  I_1 parameter index of 2nd parameter
#' @param  I_2 parameter index of 3rd parameter
#' @param  I_3 parameter index of 4th parameter
#' @param  z dimension of the parameter subspace 
#' @param  Verbose more/less (1/0) output information 
#' @param  Minimization TRUE/FALSE. 
#' @keywords mle, MacKenzie, colonization, extinction
#' @useDynLib islandPLUS
#' @export
#' @examples
#' Data1 <- lakshadweep[[1]]
#' Name_of_Factors <- c("Species","Atoll","Guild")
#' Factors <- Filter(is.factor, Data1)
#' No_of_Factors <- length(Factors[1,])
#' n <- No_of_Factors + 1
#' D1 <- as.matrix(Data1[1:nrow(Data1),n:ncol(Data1)])
#' Time <- as.double(D1[1,])
#' P1 <- as.matrix(D1[2:nrow(D1),1:ncol(D1)])
## Calculation true time vector without repetitions and counting the number of
## replicates per sampling time to build the \code{Transects} vector (via calling
## native C code from dynamic library islandPLUS.so)
#' RES <- da_Counting_Replicates_per_Time( Time )
#' Time_Vector <- RES[[1]]
#' Transects   <- RES[[2]]
#' R1 <- da_MacKenzie_mle(P1, Time_Vector, Transects,
#'                        Colonization=0.5, Extinction=0.5, Detectability=0.5,
#'                        Phi_Time_0=0.5,
#'                        Tol=1.0e-8, Verbose = 1)

da_MacKenzie_mle <- function( Data, Time, Transects, 
                             Colonization = 0.1, Extinction   = 0.1,
                             Detectability = 0.99, Phi_Time_0 = 0.5,
                             Tol = 1.0E-06, MIT = 100, 
                             C_MAX = 2.0, C_min = 0.0,
                             E_MAX = 2.0, E_min = 0.0,
                             D_MAX = 0.999, D_min = 0.001, 
                             P_MAX = 0.999, P_min = 0.001, 
                             I_0 = 0, I_1 = 1, I_2 = 2, I_3 = 3, z = 4, 
                             Verbose = 0, Minimization = TRUE )
{
    No_C = 100;
    No_E = 100;
    No_D = 100;
    No_P = 100; 
     
        S <- nrow(Data);    # No of species 
        N <- ncol(Data);    # Total No of Temporal Observations
        n <- length(Time);  # No of sampling Times, 
        Z <- 4;             # Maximum dimension of parameter space

        if (z > Z)  return("Error: Dimension of the parameter space is greater than
 the maximum allowed parameter space dimension (4)");
        
# Input argument Data should be transposed because the internal storage of
# matrices in R is just the oposite as the criterion I used when coding the shared
# object this R function calls:
	      P <- array( t(Data), c(S*N) );
        
        if( length(Time) != length(Transects) ) return("Error: Number of sampling times does not match the length of the Transects vector");
        if (sum(Transects) != N ) return("Error: total number of transects does not match the number of columns of input data matrix");
        
        C_Range <- double(2); C_Range <- c(C_min, C_MAX);
        E_Range <- double(2); E_Range <- c(E_min, E_MAX);
        D_Range <- double(2); D_Range <- c(D_min, D_MAX);
        P_Range <- double(2); P_Range <- c(P_min, P_MAX);

        Index <- integer(4);             Index <- c(I_0, I_1, I_2, I_3);
        Discretization <- integer(4);    Discretization <- c(No_C, No_E, No_P, No_P);
        Tr <- integer( n );              Tr <- Transects;

        if( Minimization ) { MINIMIZATION = 1; }
        else {               MINIMIZATION = 0; }
        
        Value <-  0; 
        res <- .C("R_SHLIB___mle_MacKenzie_NLLikelihood_Minimization", PACKAGE="islandPLUS",
                  as.double(P), as.integer(S), as.integer(N),
                  as.double(Time), as.integer(Tr), as.integer(n),
                  as.double(Colonization), as.double(C_Range), 
                  as.double(Extinction), as.double(E_Range),
                  as.double(Detectability), as.double(D_Range),
                  as.double(Phi_Time_0), as.double(P_Range),
                  as.integer(z), as.integer(Z),
                  as.integer(Index), as.integer(Discretization),
                  as.double(Tol), as.integer(MIT),
                  as.double(Value),
                  as.integer(Verbose), as.integer(MINIMIZATION) );
         
        RES <- list(C=res[[7]], E=res[[9]],
                    D=res[[11]], P=res[[13]],
                    NLL=res[[21]]);

        RES;
    }
                  
        
                       
                       
                       
