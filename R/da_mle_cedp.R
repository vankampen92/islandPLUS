#' Likelihood approach for esitmating colonization/extinction with perfect or imperfect detectability   
#' 
#' \code{da_mle_cedp} conducts maximum likelihood estimation of colonization/extinction parameters
#' of different data sets. This function can handle imperfect detectability and missing data
#' defining a heterogeneous sampling structure across input data matrix rows. The function generates, 
#' as an output, either a 3-column matrix (Colonization, Extinction, Negative LogLikelihood) or
#' 5-column matrix (Colonization, Extinction, Detectability, P_0, Negative LogLikelihood), depending
#' on the value of the input parameter PerfectDetectability (either \code{TRUE} or \code{FALSE}). 
#'
## Details:
#' The input is a data frame containing presence data per time (in cols).  Different factors (for instance,
#' OTU, location, etc) can slide the initial data frame accordingly. Model parameters will be estimated
#' for each of these groups independently that correspond to each level of the chosen factor. If Minimization
#' is 0, then no maximum likelihood estimation is performed and only the likelihood evaluation at the input model
#' parameter values is returned. Searches are based on the Nelder-Mead simplex method, but conducted in a
#' bounded parameter space which means that in case a neg loglikelihood (NLL) evaluation is called out from these
#' boundaries, the returned value for this NLL evaluation is artifically given as the maximum number the machine can
#' hold. Each group is named by a short-length-character label (ideally, 3 or 4 characters). All labels should
#' have the same character length to fulfill memmory alignment requirements of the shared object called by
#' .C(...) function. I_0, I_1, I_2, I_3 are model parameter keys. They are used to define a 4D-vector (Index). The
#' search will take place on the full parameter space defined by model parameters (I_0, I_1) if PefectDetectability
#' is \code{TRUE} or, alternatively, defined by (I_0, I_1, I_2, I_3) if PerfectDetectability is \code{FALSE}.
#' Model parameter keys correspond to colonization (0), extinction (1), detectability (2), and P_0 (3) model
#' parameters. For instance, if (I_0, I_1) is (1, 0), the search will take place whitin the paremeter space defined
#' by extinction, as the first axis, and colonization, as the second. 
## Input arguments:
#' @param  Data data frame containing presence data per time (in cols) and sites (in rows) 
#' @param  Time an array of length n containing sampling times  
#' @param  Tags array of names (one short for each level of the factor analyzed)  
#' @param  Factor column number containing the 'data frame' factor used to split total data into level-based groups
#' @param  Colonization guess value to initiate search / parameter value
#' @param  Extinction guess value to initiate search / parameter value
#' @param  Detectability_Value guess value to initiate search / parameter value
#' @param  Phi_Time_0_Value guess value to initiate search / parameter value
#' @param  Tol stopping criteria of the search algorithm. 
#' @param  MIT max number of iterations of the search algorithm.
#' @param  C_min min value for the possible range of colonization values
#' @param  C_MAX max value for the possible range of colonization values
#' @param  E_min min value for the possible range of colonization values
#' @param  E_MAX max value for the possible range of colonization values
#' @param  D_min min value for the possible range of colonization values
#' @param  D_MAX max value for the possible range of colonization values
#' @param  P_min min value for the possible range of colonization values
#' @param  P_MAX max value for the possible range of colonization values   
#' @param  I_0 has to be 0 or 1. Defaults to 0 
#' @param  I_1 has to be 0 or 1. Defaults to 1
#' @param  I_2 has to be 0, 1, 2, or 3. Defaults to 2 
#' @param  I_3 has to be 0, 1, 2, or 3. Defaults to 3
#' @param  Minimization, 1/0. If Minimization is 0, then no minimization is performed. 
#' @param  z dimension of the parameter subspace for which the optimization process will take place. Defaults to 2 
#' @param  Verbose more/less (1/0) information about the optimization algorithm will be printed out.
#' @param  MV_FLAG missing Value Flag (to specify sites and times where no sample exists)
#' @param  PerfectDetectability TRUE means 'Perfect Detectability'. Of course, FALSE means 'Imperfect Detectability' 
#' @keywords maximum likelihood, colonization, extinction
#' @useDynLib islandPLUS
#' @export     
#' @examples
## Example (Lakshadweep Islands --- Imperfect Detecability)
#' Data <- lakshadweepPLUS[[1]]
#' Guild_Tag = c("Alg","Cor","Mac","Mic","Omn","Pis","Zoo")
#' Time <- as.vector(c(2000,2000,2001,2001,2001,2001,2002,2002,2002,2002,2003,2003,2003,2003,2010,2010,2011,2011,2011,2011,2012,2012,2012,2012,2013,2013,2013,2013))
#' R <- da_mle_cedp(Data, Time, Factor=3, Tags=Guild_Tag, PerfectDetectability=FALSE, z=4)
#' Guild_Tag = c("Agt","Kad","Kvt")
#' R <- da_mle_cedp(Data, Time, Factor=2, Tags=Guild_Tag, PerfectDetectability=FALSE, z=4)

da_mle_cedp <- function( Data, Time, Factor, Tags,
                       Colonization = 1.0, Extinction = 1.0,
                       Detectability_Value = 0.5, Phi_Time_0_Value = 0.5, 
                       Tol = 1.0E-08, MIT = 100, 
                       C_MAX = 10.0, C_min = 0.0,
                       E_MAX = 10.0, E_min = 0.0, 
                       D_MAX = 0.99, D_min = 0.0,
                       P_MAX = 0.99, P_min = 0.01,
                       I_0 = 0, I_1 = 1, I_2 = 2, I_3 = 3,
                       z = 2,
                       Minimization = 1, 
                       Verbose = 0, MV_FLAG = 0.1,
                       PerfectDetectability = TRUE )
{
    No_C = 100;
    No_E = 100;
    No_D = 100;
    No_P = 100;
    
    MINIMIZATION = Minimization; 
  # S: Number of OTUs (levels in Factor) considered in this classificion. 
  #    OUTs can be species, guilds, phylogentic groups, etc 
    S <- length(Tags)
    L <- split(Data, Data[[Factor]])
    Sites <- integer(S)
    No_of_Times <- integer(S);
    Nc <- length(Time);       # Nc, No of sampling Times,
    N = 0;
    for (i in 1:S) 
    {
        Sites[i] <- nrow(L[[i]]);
        No_of_Times[i] = Nc;
        N <- N + Sites[i];
    }
    if( N != nrow(Data) )
        return("Error: Total number of rows in data frame Data differs from the output of nrow(Data)");
  # N is now the total number of rows across all data matrices 
  
  # Notice here how the array function fills Time_Vector
  # with the elemens of Time in a cyclic manner until 
  # reaching a length given by its second argument. 
    Time_Vector  <- array( Time, S*Nc );
  
    Factors <- Filter(is.factor, Data);
    No_of_Factors <- length(Factors[1,]);
    n <- No_of_Factors + 1
  
  # Merging Data frames in the right order and converting the 
  # resulting data frame into a matrix also in the right ordering. 
  # Notice that the merge command should take both sort=FALSE 
  # and all=TRUE options!!!
    Data_M <- merge(L[[1]],L[[2]],all=TRUE,sort=FALSE);
    if( S > 2 ) for (i in 3:S)  Data_M <- merge(Data_M, L[[i]],all=TRUE,sort=FALSE);
    D <- as.matrix(Data_M[1:nrow(Data_M),n:ncol(Data_M)]);
    P <- array( t(D), N*Nc );
    if (PerfectDetectability == TRUE) {  
        Z <- 2;     # Maximum dimension of parameter space (Colonization, Extinction)
        if (z != Z) return("Error: Dimension of the model parameter space should be 2,  
                      that is, the colonization-extinction 2D plane.");
  
        C_Range <- double(2); C_Range <- c(C_min, C_MAX);
        E_Range <- double(2); E_Range <- c(E_min, E_MAX);
  
        Index <- integer(4);             Index <- c(I_0, I_1, 2, 3);
        Discretization <- integer(4);    Discretization <- c(No_C, No_E, 1, 1);
  
        Result <- double(3*S);
        Res <- .C("R_SHLIB___mle_NLLikelihood_Minimization", # PACKAGE="islandPLUS",
                  as.integer(S), as.character(Tags),
                  as.double(P), as.integer(Sites), 
                  as.double(Time_Vector), as.integer(No_of_Times),
                  as.double(MV_FLAG),
                  as.double(Colonization), as.double(C_Range), 
                  as.double(Extinction), as.double(E_Range),
                  as.integer(z), as.integer(Z),
                  as.integer(Index), as.integer(Discretization),
                  as.double(Tol), as.integer(MIT),
                  as.integer(Verbose), as.integer(MINIMIZATION),
                  as.double(Result));
  
  # Result corresponds to the 20th argument of the .C(...) call list. 
  # True arguments are counted from 1 to the last. The first two args are not
  # considered in this count.     
  # Arg 0 is R_SHLIB___mle_NLLikelihood_Minimization, this is, the name 
  # of the shared object.
        Result  <- Res[[20]];
        
        R <- t(array(Result, c(3, S)));
    }
    else {
        Z <- 4;     # Maximum dimension of parameter space (Col, Ext, Dtc, P_0)
        if (z != Z) return("Error: Dimension of the model parameter space should be 4,  
                      that is, the colonization-extinction-detectability-P_0 hypervolum.");

        C_Range <- double(2); C_Range <- c(C_min, C_MAX);
        E_Range <- double(2); E_Range <- c(E_min, E_MAX);
        D_Range <- double(2); D_Range <- c(D_min, D_MAX);
        P_Range <- double(2); P_Range <- c(P_min, P_MAX);
        
        Index <- integer(4);             Index <- c(I_0, I_1, I_2, I_3);
        Discretization <- integer(4);    Discretization <- c(No_C, No_E, No_D, No_P);
        
        Result <- double(5*S);
        Res <- .C("R_SHLIB___mle_MacKenzie_NLL_Uneven_Minimization",
                  as.integer(S), as.character(Tags),
                  as.double(P), as.integer(Sites), 
                  as.double(Time_Vector), as.integer(No_of_Times),
                  as.double(MV_FLAG),
                  as.double(Colonization), as.double(C_Range), 
                  as.double(Extinction), as.double(E_Range),
                  as.double(Detectability_Value), as.double(D_Range), 
                  as.double(Phi_Time_0_Value), as.double(P_Range),
                  as.integer(z), as.integer(Z),
                  as.integer(Index), as.integer(Discretization),
                  as.double(Tol), as.integer(MIT),
                  as.integer(Verbose), as.integer(MINIMIZATION),
                  as.double(Result));

  # Result corresponds to the 24th argument of the .C(...) call list. 
  # True arguments are counted from 1 to the last. The first two args are not
  # considered 'true args' in this count.         
  # Arg 0 is MODEL_SELECTION_UPGMA_R_FUNCTION, this is, the name 
  # of the shared object.
        Result  <- Res[[24]];
        
        R <- t(array(Result, c(5, S)));
    }

    R
}
                  
        
                       
                       
                       
