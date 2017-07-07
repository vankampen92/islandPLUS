#include "HEADERS.h"
#include "da_IBT_Functions.h"

void Rates_into_Probabilities( double c, 
			       double e,
			       double * T10,
			       double * T01, 
			       double DT )
{
  /* Input:
     . e, Extinction_Rate
     . c, Colonization_Rate
     . DT, Time interval
     
     Output:
     . T10, Colonization probability
     . T01, Extinction probability
  */

  (*T10) = c/(e+c)* ( 1.0 - exp(-(e+c)*DT) );
  (*T01) = e/(e+c)* ( 1.0 - exp(-(e+c)*DT) );
}
