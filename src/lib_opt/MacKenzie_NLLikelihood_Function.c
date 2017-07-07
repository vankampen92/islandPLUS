#include <MODEL.h>

void program_exit_for_matrix_corruption(int Sp, int Presence, int Time)
{
  Rprintf("Species: %d\t State: %d\t Time: %d\n", Sp, Presence, Time); 
  Rprintf("State can only be either 0 or 1\n");
  Rprintf("The program will exit\n");
  error(0,0,"Program aborted");
}

int Determining_Total_No_of_Absences (double * Presence_Data, int N, 
				      int * Transects, int * Time_Index, int n)
{
  int No_of_Observed_Absences;
  int i,j,m, n_0,k;
      
  n_0=0; k = 0;
  for(i=0; i<n; i++) {
    m = 0;
    for (j=n_0; j<(n_0+Transects[i]); j++) {
      if ( j >= N ) Rprintf(" j = %d\t N = %d\n", j, N); if( j >= N ) error(0,0,"Program aborted");
      if (Presence_Data[j] == 0.0) m++;
    }
    if( m == Transects[i] ) Time_Index[k++] = i; 
    
    n_0 += Transects[i];
  }
  
  No_of_Observed_Absences = k;
  return(No_of_Observed_Absences);
}

void True_Presence_Absence_Data( int * Presence, int T,
				 int * Binary, int * Time_Index, 
				 int No_of_Absences )
{
  int i;
  
  for( i = 0; i<T; i++ ) Presence[i] = 1;
  for( i = 0; i<No_of_Absences; i++ ) Presence[Time_Index[i]] = Binary[i];  
}

double MacKenzie_NLLikelihood_Calculation ( double ** Presence_Data, int No_of_SPECIES, int N, 
					    double * Time_Vector, int * Transects, int T, 
					    double Colonization_Rate, 
					    double Extinction_Rate,
					    double Detectability_Value, 
					    double Phi_0 )
{
  /* Input arguments:
     . Presence_Data[][] : 1/0 Data matrix. 
     . No_of_SPECIES: Number of rows of Data matrix. 
     . N: Number of columns of the presence data matrix.
     . Time_Vector : array containing sampling times. 
     . Transects: array containing number of transects per sampling time.
     . T : Number of empirical observations (sampling times)
     . Colonization_Rate: colonization rate in T^{-1} units
     . Extintion_Rate: extinction rate in T^{-1} units
     . Detectability_Value: detection probability. 
     . Phi_0: Probability of species presence at time 0

     Output arguments:
     . Neg Log Likelihood (see MacKenzie 2003)
  */ 
  int i,j,k,m,n,n_0;
  int k_0, k_1;
  long double X, X_TOTAL;
  long double Y;
  double NLL; 
  int No_of_Absences;
  int No_of_Compatible_Dynamics; 

  double ** TM = (double **)calloc(2, sizeof(double *) );
  TM[0]= (double *)calloc(2, sizeof(double) );
  TM[1]= (double *)calloc(2, sizeof(double) );

  int PAR_DEFINITION; 
  /* Checking if parameters are well defined. If not 
     the likelihood makes no sense and is not calculted. 
     The return value is only a very large number */
  PAR_DEFINITION = Checking_for_Parameter_Correctness( Colonization_Rate, 
						       Extinction_Rate,
						       Detectability_Value, 
						       Phi_0 );	
					
  if ( PAR_DEFINITION  ==  0 ) {
   
    NLL = DBL_MAX;
   
    return (NLL);
  } 
  else {
    /* B E G I N : Main likelihood evaluation starts */
    X_TOTAL = 0.0; 
    for(j=0; j<No_of_SPECIES; j++) {
  
      int * Time_Index   = (int *)calloc(T, sizeof(int) );
      No_of_Absences = Determining_Total_No_of_Absences (Presence_Data[j], N,
							 Transects, Time_Index, T); 
      
      No_of_Compatible_Dynamics = (int)pow( 2.0, (double)No_of_Absences );
      
      int ** Binary_Combination = (int **)calloc( No_of_Compatible_Dynamics, 
						  sizeof(int *) );
      for (k=0; k<No_of_Compatible_Dynamics; k++) 
	Binary_Combination[k] = (int *)calloc(No_of_Absences, sizeof( int ) );
      
      Create_Binary_Combination( Binary_Combination, No_of_Compatible_Dynamics, 
			       No_of_Absences );
  
      X = 0.0;
      for(k=0; k<No_of_Compatible_Dynamics; k++) {
	
	int * Presence     = (int *)calloc(T, sizeof(int) );
	True_Presence_Absence_Data( Presence, T,
				    Binary_Combination[k], Time_Index, 
				    No_of_Absences );
	Y = 0.0; 
	if (Presence[0] == 0) {
	  m = 0;
	  for(n=0; n<Transects[0]; n++) if (Presence_Data[j][n] == 0) m++;
	  if (m != Transects[0]) error(0,0,"Number of Transects does not match");
	
	  Y += logl( 1.0 - (long double)Phi_0 ); 
	}
	else if (Presence[0] == 1) {
	  m = 0;
	  for(n=0; n<Transects[0]; n++) if (Presence_Data[j][n] == 0) m++;
	  
	  Y += (long double)m * logl( 1.0 - (long double)Detectability_Value );
	  
	  Y += (long double)(Transects[0]-m) * logl( (long double)Detectability_Value );
	 
	  //Y += log(Combinarory(m, Transects[0]));
	
	  Y += logl ( (long double)Phi_0 );
	}
	else {
	  program_exit_for_matrix_corruption(j, Presence[0], 0);
	}
	
	n_0 = Transects[0];
	for(i=1; i<T; i++) {
	
	  double Delta_T = Time_Vector[i] - Time_Vector[i-1];
	  Transition_Matrix( TM, 2,2,
			     Colonization_Rate,
			     Extinction_Rate,
			     Delta_T );
	  
	  /* Transiction: k_0 ---> k_1 */
	  k_0 = Presence[i-1];    if( k_0 != 0 && k_0 != 1 ) error(0,0,"Program aborted");  
	  k_1 = Presence[i];      if( k_1 != 0 && k_1 != 1 ) error(0,0,"Program aborted");  
	
	  Y += logl( (long double)TM[k_1][k_0] ) ;
	  
	  if( Presence[i] == 0 ) {
	    Y += 0.0; // log(1.0) = 0.0;
	    m = 0;
	    for(n = n_0; n<(n_0+Transects[i]); n++) if (Presence_Data[j][n] == 0) m++;
	    if ( m != Transects[i] ) error(0,0,"Number of Transects do not match");
	  }
	  else if (Presence[i] == 1) {
	    m = 0;
	    for(n = n_0; n<(n_0+Transects[i]); n++) if (Presence_Data[j][n] == 0) m++;
	    
	    Y += (long double)m * logl(1.0 - (long double)Detectability_Value );
	    
	    Y += (long double)(Transects[i]-m) * logl( (long double)Detectability_Value ) ;

	    if(Transects[i] < m) error(0,0,"Program aborted");
	    //X += log(Combinarory(m, Transects[0]));
	  }
	  else {
	    program_exit_for_matrix_corruption(j, Presence[i], i);
	  }
	  
	  n_0 += Transects[i];
	  if( i == (T-1) ) if( n_0 != N ) error(0,0,"Program aborted");
	}
	free(Presence);
	
	/* If Y is too high, X should take the value of 0.0 */
	X += expl(Y);

	if( Y >= 0.0 ) error(0,0,"Program aborted");
      }

      free(Time_Index);
      for (k=0; k<No_of_Compatible_Dynamics; k++) free(Binary_Combination[k]);
      free(Binary_Combination);
      
      if ( X == 0.0 )      X_TOTAL += DBL_MAX/(double)No_of_SPECIES;
      else                 X_TOTAL += -logl(X);
    }
    
    free(TM[0]); free(TM[1]);
    free(TM);
    
    NLL = (double)X_TOTAL; 
    return(NLL);
    /*     E N D : Main likelihood evaluation ends!!! */
  }
}
 
	

