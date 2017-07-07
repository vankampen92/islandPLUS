#' Auxiliary function for MacKenzie mle methods  
#'
#' \code{da_Counting_Replicates_per_Time} counts transects/replicates per sampling time.
#' The number of replicates per sampling time can perfectly differ.
#' 
# Input arguments:
#' @param Time an array containing sampling times (including time repetitions)
#'
#' The output of this function generates a list containing a time vector 'Time_Vector'
#' with all sampling times (without repetitions) and an integer vector 'Transects'
#' counting the number of repetitions or replicates per sampling time.
#' @keywords MacKenzie, colonization, extinction, detectability
#' @useDynLib islandPLUS
#' @export     
#' @examples
#' Time <- c(2000,2000,2000,2001,2001,2002,2002,2002,2002)
#' R <- da_Counting_Replicates_per_Time( Time )

da_Counting_Replicates_per_Time <- function( Time ) 
{
  N <- length(Time);       # No of sampling Times, 
  
  Time_Vector <- double(N);

  Transects <- integer(N);
  
  No_of_Times <- N; 
  
  res <- .C("Counting_Replicates_per_Time", PACKAGE="islandPLUS",
            as.double(Time), as.integer(N),
            as.double(Time_Vector), as.integer(Transects),   
            as.integer(No_of_Times) );
  
  # Result corresponds to the 20th argument of the .C(...) call list. 
  # Arguments are counted from 0 to the last. 
  # Arg 0 is Counting_Replicates_per_Time, this is, the name 
  # of the shared object.

  Time_Vector <- double(res[[5]]);
  Transects <- integer(res[[5]]);
  
  for( i in 1:res[[5]] ) {
      Time_Vector[i] <- res[[3]][i]
      Transects[i] <- res[[4]][i]
  }
  
  Result  <- list(Time_Vector, Transects);

  Result; 
}
                  
        
                       
                       
                       
