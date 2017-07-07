#include "HEADERS.h"
#include "da_IBT_Functions.h"

// This two functions below belong to standard I/O procedures. 
// (see IO_Procedures_Standard.h for function prototypes)
// They provide general procedures to filtering out missing values, 
// when they are initially marked with a given flag. 

void IO_Filtering_Out_Missing_Values ( int No_of_SPECIES,  
				       double *** Presence, int * No_of_SITES, 
				       double ** Time_Vector, int * No_of_TIMES, 
				       double *** Sp_Time, int ** No_of_Sp_Times, 
				       double MISSING_VALUE_FLAG )
{
  /* This function is only a wrapper to IO_Filtering_Out_Matrix (...). Essentially, 
     it applies this function over a loop that goes over each intial Presence matrix.  
  */
  int i; 
  int ACTUAL_No_of_ROWS;
  
  for(i=0; i<No_of_SPECIES; i++){
    ACTUAL_No_of_ROWS = No_of_SITES[i];
 
    IO_Filtering_Out_Matrix( Presence[i], &ACTUAL_No_of_ROWS, 
			      Time_Vector[i], No_of_TIMES[i],    
			      Sp_Time[i], No_of_Sp_Times[i], 
			      MISSING_VALUE_FLAG );
    
    No_of_SITES[i] = ACTUAL_No_of_ROWS;  
  }
}

void IO_Filtering_Out_Matrix( double ** Presence, int * No_of_SPECIES,  
			      double * Time, int No_of_TIMES,    
			      double ** Sp_Time, int * No_of_Sp_Times, 
			      double MISSING_VALUE_FLAG )
{
  /* This function calculates true time vectors and their true length.
     from a matrix Presence, which has a flag value (MISSING_VALUE_FLAG) 
     to mark those times that are missing. 

     The output is a Sp_Time array that gives the specific time vector 
     (Sp_Time) and its actual length (No_of_Sp_Times) for row of the 
     initial matrix once missing values have been discounted. 

     The function also counts the number of actual rows of the initial 
     Presence matrix. 

     Input arguments:
     . Presence[][]: presence/absence values (including missing values)
     . (* No_of_SPECIES) total number of rows in the input Presence 
     matrix. 
     . Time[]: sampling times corresponding to each column of the initial
     Presence matrix
     . No_of_TIMES: total number of columns of the input Presence matrix.
     . MISSING_VALUE_FLAG: a double to mark missing values 

     Output arguments:
     . Sp_Time[i][j]: j-th sampling time associated to the i-th row 
     once missing values in Presence has been excluded. 
     . No_of_Sp_Times[i]; number of sampling times associated to the 
     i-th row once missing values in Presence has been excluded. 
     . (* No_of_SPECIES) total number of effective rows in the input 
     Presence matrix for which there is, at least, two times with 
     non-missing sampling times. 
     . Presence[][] the output is the same as the initial matrix once 
     deffective rows (those with only one sampled time) and missing values 
     have been filtered out. 
  */

  int i,j;
  int m,n;
  int S; 
  int JUMPS;

  S = (* No_of_SPECIES); 
  m = 0;
  for( i=0; i<S; i++) {
    
    JUMPS = 0; 
    for( j=0; j<No_of_TIMES; j++ ) 
      if( Presence[i][j] != MISSING_VALUE_FLAG ) JUMPS++;
    
    n = 0;
    if( JUMPS > 1 ) { 
      for( j=0; j<No_of_TIMES; j++ ) { 
	if( Presence[i][j] != MISSING_VALUE_FLAG ) {
	  Presence[m][n] = Presence[i][j];
	  Sp_Time[m][n]  = Time[j];
	  
	  if (Presence[m][n] != 0.0 && Presence[m][n] != 1.0 ) error(0,0,"Program aborted");
	  
	  n++;
	}
      }
      No_of_Sp_Times[m] = n;
      if (n > 0) m++;
    }
  }
  
  (* No_of_SPECIES) = m;
}
