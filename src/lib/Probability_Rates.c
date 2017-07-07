#include "HEADERS.h"
#include "da_IBT_Functions.h"

void Probability_Rates( double C, 
			double E, 
			double * c, 
			double * e,
			double DT )
{
  /* Input:
     . E, Extinction_Probability in  discrete time intervals DT
     . C, Colonization_Probability in  discrete time intervals DT
     . DT, Time interval
     
     Output:
     . e, Extinction_Rate (instantaneous rate)
     . c, Colonization_Rate (instantaneous rate)
  */
  
  double R_E, R_C;

  R_E = E/C;
  R_C = C/E;
  
  // double l = log( 1.0 - (C+E) );
  // printf (" puto log( 1.0 - (C+E) ) = %g\n", l );
  // getchar();

  (* c) = - 1.0/DT * log(1.0 - (C+E)) /( 1.0 + R_E ); 
  (* e) = - 1.0/DT * log(1.0 - (C+E)) /( 1.0 + R_C ); 
}

