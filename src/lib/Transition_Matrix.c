#include "HEADERS.h"
#include "da_IBT_Functions.h"

void Transition_Matrix( double **  T, int Nr, int Nc, 
			       double c, 
			       double e,
			       double DT )
{
  double T00, T01, T10, T11; 
  /* Input:
     . e, Extinction_Rate
     . c, Colonization_Rate
     . Nr, # Rows. 
     . Nc, # Columns. 
     . DT, Time interval
     
     Output:
     . T. Transition Matrix:
  */
  T[0][0] = T00 = 1.0 - c/(e+c)*( 1.0 -exp(-(e+c)*DT) );
  T[1][0] = T10 = c/(e+c)*(1.0 -exp(-(e+c)*DT) );
  T[0][1] = T01 = e/(e+c)*(1.0 -exp(-(e+c)*DT));
  T[1][1] = T11 = 1.0 - e/(e+c)*(1.0 -exp(-(e+c)*DT));

  // printf("T00 = %g, T10 = %g, T01 = %g, T11 = %g\n", T00, T10, T01, T11);
}


