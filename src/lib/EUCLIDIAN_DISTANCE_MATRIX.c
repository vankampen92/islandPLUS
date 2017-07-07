#include "HEADERS.h"
#include "da_IBT_Functions.h"

void EUCLIDIAN_DISTANCE_MATRIX( double * C, double * E, int N,  
				double ** D )
{
  /* This function calculates the Euclidian distance between pairs on 
     the colonization-extinction parameters space.
  */
  int i,j;
  
  for(i=1; i<N; i++) 
    for(j=0; j<i; j++) 
      D[i][j] = sqrt( (C[i] - C[j])*(C[i] - C[j]) + (E[i] - E[j])*(E[i] - E[j]) );
}

void EUCLIDIAN_DISTANCE_MATRIX_MacKENZIE ( double * C, double * E,
					   double * D, double * P,
					   int N, 
					   double ** Distance )
{
  /* This function calculates the Euclidian distance between pairs on 
     the colonization-extinction parameters space.
  */
  int i,j;
  
  for(i=1; i<N; i++) 
    for(j=0; j<i; j++) 
      Distance[i][j] = sqrt( (C[i] - C[j])*(C[i] - C[j]) + (E[i] - E[j])*(E[i] - E[j]) + (D[i] - D[j])*(D[i] - D[j]) + (P[i] - P[j])*(P[i] - P[j]));

}
