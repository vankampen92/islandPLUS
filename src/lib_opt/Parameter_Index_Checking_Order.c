#include <MODEL.h>

void Parameter_Index_Checking_Ordering (int * A, int * D, int N, 
					int * No_C, int * No_E, int * No_D, int * No_P)
{
  /* Input parameters: 
     . A[i]: Definition of the parameter space through a list of integers
     containing the keys corresponding to model parameters.  
     In general, the array A should be a permuation of the first N integers. 
     If this is not true, the program aborts. 
     . N: Length of arrays A[] and D[]: Maximum number of model parameters. 
     . D: Parameter Space Discretization: No of discrete points in each direction of
          the suitably redefined (sub)parameter space. 
	  
     Output parameters: 
     . No_C: Discretization in the first (key = 0) parameter
     . No_E: Discretization in the second (key = 1) parameter
     . No_D: Discretization in the third (key = 2) parameter
     . No_P: Discretization in the fourth (key = 3) parameter
  */   
  int i,j,n,m;
  double * Index = (double *)calloc( N, sizeof(double) );
  size_t    * p     = (size_t   *)calloc( N, sizeof(size_t)    );
  
  m = 0;
  for(i=0; i<N; i++){
    n = 0;
    for (j=0; j<N; j++) if(A[j] == i) n++; 
    if (n == 1) m++;

    Index[i] = (double)A[i];
  }
  if( m != N ) error(0,0, "Program aborted");
    
  /* gsl_sort_index (...):
       This function indirectly sorts the n elements of the array 'Index' 
       with stride 'double' into ascending order, storing the resulting permutation 
       in p. The array p must be allocated with a sufficient length to store the 
       n elements of the permutation. The elements of p give the index of the array 
       element which would have been stored in that position if the array had been 
       sorted in place. The array Index is not changed. 
  */
    gsl_sort_index (p, Index, sizeof(double), N);
    
    (* No_C) = D[p[0]]; 
    (* No_E) = D[p[1]]; 
    (* No_P) = D[p[2]];
    (* No_D) = D[p[3]];
    
    free(Index); free(p);
}
