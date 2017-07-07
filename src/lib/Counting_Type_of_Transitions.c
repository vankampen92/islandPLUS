#include "HEADERS.h"
#include "da_IBT_Functions.h"
     
void Counting_Type_of_Transitions( double ** Presence, int S, 
				   double ** Time, int * T, 
				   int * N00, int * N01, int * N10, int * N11 )
 {
      int i,j;
      int n00, n01, n10, n11; 
      
      n00 = n01 = n10 = n11 = 0;
      for( i=0; i<S; i++ ) {
	for(j=1; j<T[i]; j++ ) {
	  if (     Presence[i][j-1] == 0.0 && Presence[i][j] == 0.0) { n00++; }
	  else if (Presence[i][j-1] == 0.0 && Presence[i][j] == 1.0) { n10++; }
	  else if (Presence[i][j-1] == 1.0 && Presence[i][j] == 0.0) { n01++; }
	  else if (Presence[i][j-1] == 1.0 && Presence[i][j] == 1.0) { n11++; }
	  else    { Rprintf(" Error in Counting_Type_of_Transitions(...)\n"); error(0,0,"Program aborted"); }
	}
      }
      
      (* N00 ) = n00; 
      (* N01 ) = n01;
      (* N10 ) = n10; 
      (* N11 ) = n11; 
 }



