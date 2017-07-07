void EUCLIDIAN_DISTANCE_MATRIX( double * C, double * E, int N,  
				double ** D );

void EUCLIDIAN_DISTANCE_MATRIX_MacKENZIE ( double * C, double * E,
                                           double * D, double * P,
                                           int N,
                                           double ** Distance ); 

void Probability_Rates( double C,
                        double E,
                        double * c,
                        double * e,
                        double DT );

void Rates_into_Probabilities( double c, 
			       double e,
			       double * T10,
			       double * T01, 
			       double DT );

void Transition_Matrix( double **  T, int Nr, int Nc,
			double c,
			double e,
			double DT );


void int_buffer_rec(int ** Number_List, int N, 
		    int * number, int n, int length); 

void Create_Binary_Combination( int ** Binary_Combination, int N, 
				int LENGTH ); 

void Parameter_Index_Checking_Ordering (int * A, int * D, int N, 
					int * No_C, int * No_E, int * No_D, int * No_P); 

#include  "IO_Missing_Values.h"

#include "IO_AKAIKE_Model_Selection.h"

/* Counting type of transiton fo heterogeneous input data matrix */
void Counting_Type_of_Transitions( double ** Presence, int S, 
				   double ** Time, int * T, 
				   int * N00, int * N01, int * N10, int * N11 );

/* Function to count replicates per time required to work with 
   MacKenzie functions */
void Counting_Replicates_per_Time(double * Vector, int * Temporal_Observations,
                                  double * Time_Vector, int * Transects,
                                  int * No_of_TIMES);

/*     E N D : --------------------------------------------------*/

/* B E G I N : Prototype definition UPGMA algorithm -----------*/
#include "upgma_clustering.h"
/*     E N D : ------------------------------------------------*/
