#include <MODEL.h>

void MODEL_SELECTION_UPGMA_DRIVER ( int No_of_SPECIES, char ** Species_Tag, 
				    double *** Presence, int * No_of_SITES, 
				    double ** Time_Vector, int * No_of_TIMES, 
				    double MISSING_VALUE_FLAG, 
				    double Colonization_Rate, double * C_Range, 
				    double Extinction_Rate, double * E_Range, 
				    int * No_of_PARAMETERS,
				    int * No_of_PARAMETERS_MAX,
				    int * Index, 
				    int * Discretization,
				    double * Tolerance, 
				    int * No_of_ITERATIONS,
				    int * Verbose, 
				    int * Minimization,
				    double ** Results );

void MODEL_SELECTION_UPGMA_R_FUNCTION( int * pNo_of_SPECIES, char ** Species_Tag, 
				       double * R_Presence, int * No_of_SITES, 
				       double * R_Time_Vector, int * No_of_TIMES, 
				       double * MISSING_VALUE_FLAG, 
				       double * Colonization_Rate, double * C_Range,
				       double * Extinction_Rate, double * E_Range, 
				       int * No_of_PARAMETERS,
				       int * No_of_PARAMETERS_MAX,
				       int * Index, 
				       int * Discretization,
				       double * Tolerance, 
				       int * No_of_ITERATIONS,
				       int * Verbose, 
				       int * Minimization,
				       double * R_Results )
{ 
  int i,j,k,l,s,  n;
  int No_of_SPECIES; 
  No_of_SPECIES = (* pNo_of_SPECIES);

  double ** Results = (double **)calloc( No_of_SPECIES, sizeof(double *) );
  for( i=0; i<No_of_SPECIES; i++ ){
    Results[i] = (double *)calloc( 6, sizeof(double) );
  }  
  double *** Presence = (double ***)calloc(No_of_SPECIES, sizeof(double **));
  double ** Time_Vector = (double **)calloc(No_of_SPECIES, sizeof(double *));
  for( i=0; i<No_of_SPECIES; i++ ){    
    Time_Vector[i] = (double *)calloc(No_of_TIMES[i], sizeof(double));
    Presence[i] = (double **)calloc(No_of_SITES[i], sizeof(double *));
    for (j=0; j<No_of_SITES[i]; j++)
      Presence[i][j] = (double *)calloc(No_of_TIMES[i], sizeof(double ));
  }

  n=0; 
  for(i=0; i<No_of_SPECIES; i++) 
    for(j=0; j<No_of_TIMES[i]; j++) 
      Time_Vector[i][j] = R_Time_Vector[n++];
  
  int No_of_COLUMNS; 
  No_of_COLUMNS = No_of_TIMES[0]; 
  for ( i=0; i<No_of_SPECIES; i++ )   
    if ( No_of_COLUMNS != No_of_TIMES[i] ) error(0,0,"Program aborted");  
  Rprintf( "Number of Columns: %d\n", No_of_COLUMNS );

  n=0; 
  for(i=0; i<No_of_SPECIES; i++) 
    for(j=0; j<No_of_SITES[i]; j++)
      for(k=0; k<No_of_TIMES[i]; k++)
	Presence[i][j][k] = R_Presence[n++];

  MODEL_SELECTION_UPGMA_DRIVER (No_of_SPECIES, Species_Tag, 
				Presence, No_of_SITES, 
				Time_Vector, No_of_TIMES, 
				(* MISSING_VALUE_FLAG), 
				(* Colonization_Rate), C_Range,  
				(* Extinction_Rate), E_Range, 
				No_of_PARAMETERS, No_of_PARAMETERS_MAX,
				Index, Discretization, 
				Tolerance, No_of_ITERATIONS, 
				Verbose, Minimization, 
				Results); 
  n=0;
  for(i=0; i<No_of_SPECIES; i++) 
    for(j=0; j<6; j++)
      R_Results[n++] = Results[i][j];
  
  for( i=0; i<No_of_SPECIES; i++ ) free(Results[i]);
  free(Results);

  for( i=0; i<No_of_SPECIES; i++ ){
    free(Time_Vector[i]);
    for (j=0; j<No_of_SITES[i]; j++) free( Presence[i][j] );
    free(Presence[i]); 
  }
  free( Presence ); free( Time_Vector); 
  /*   END: ------------------------------------------------------------ 
   */
}

void MODEL_SELECTION_UPGMA_DRIVER ( int No_of_SPECIES, char ** Species_Tag, 
				    double *** Presence, int * No_of_SITES, 
				    double ** Time_Vector, int * No_of_TIMES, 
				    double MISSING_VALUE_FLAG, 
				    double Colonization_Rate, double * C_Range, 
				    double Extinction_Rate, double * E_Range, 
				    int * No_of_PARAMETERS,
				    int * No_of_PARAMETERS_MAX,
				    int * Index, 
				    int * Discretization,
				    double * Tolerance, 
				    int * No_of_ITERATIONS,
				    int * Verbose, 
				    int * Minimization,
				    double ** Results )
{
  int i,j,k,l,s,  n;
  
  double *** Sp_Time = (double ***)calloc(No_of_SPECIES, sizeof(double **));
  for( i=0; i<No_of_SPECIES; i++ ){
    Sp_Time[i] = (double **)calloc(No_of_SITES[i], sizeof(double *) );
    for (j=0; j<No_of_SITES[i]; j++) { 
      Sp_Time[i][j] = (double *)calloc(No_of_TIMES[i], sizeof(double ));
      for(k = 0; k<No_of_TIMES[i]; k++) {
	Sp_Time[i][j][k] = Time_Vector[i][k];
      }
    }
  }
  int ** No_of_Sp_Times = (int **)calloc(No_of_SPECIES, sizeof(int *));
  for( i=0; i<No_of_SPECIES; i++ ) {
    No_of_Sp_Times[i] = (int *)calloc(No_of_SITES[i], sizeof(int) );
    for( j=0; j<No_of_SITES[i]; j++) {
      No_of_Sp_Times[i][j] = No_of_TIMES[i];
    }
  }
  /*
    Before filtering out missing values, full data matrices and sampling 
    times are allocated in an array of length No_of_SPECIES of type 
    SP_Matrix_Data (see definition in MODEL_SP_Matrix_Data_STRUCT.h). This  
    generic data structure allows to store all these data information in 
    a convenient manner. The reason for that is that, after filetering out 
    missing values, Presence matrices are changed on place. 
  */
  SP_Matrix_Data ** Data = (SP_Matrix_Data **)calloc( No_of_SPECIES, 
						      sizeof(SP_Matrix_Data *) );
  /* Warning: No_of_COLUMNS for all data files should match!!!  */
  int No_of_COLUMNS = No_of_TIMES[0]; 
  for(i=0; i<No_of_SPECIES; i++) if( No_of_COLUMNS != No_of_TIMES[i] ) error(0,0,"Program aborted");
 
  int * Dummy = (int *)calloc(No_of_COLUMNS, sizeof(int)); //Dummy vector! 
  //No use in this context. It could handle different replicates (transects)
  //for each sampling time. 
  for( i=0; i<No_of_SPECIES; i++ ) {
    Data[i] = SP_Matrix_Data_Alloc( No_of_SITES[i], No_of_COLUMNS,
				    No_of_COLUMNS ) ;
    SP_Matrix_Data_Setup( No_of_SITES[i], No_of_COLUMNS, No_of_COLUMNS, 
			  Data[i], Presence[i], 
			  Time_Vector[i], Sp_Time[i], No_of_Sp_Times[i], 
			  Dummy, Species_Tag[i] ); 
  }
  IO_Filtering_Out_Missing_Values ( No_of_SPECIES,  
				    Presence, No_of_SITES, 
				    Time_Vector, No_of_TIMES, 
				    Sp_Time, No_of_Sp_Times, 
				    MISSING_VALUE_FLAG );
  /*       E N D :------------------------------------------------------
   */
  
  /* B E G I N : -------------------------------------------------------
     2: Calculation of mle of extinction and colonization parameters:
     a pair of parameters per each of the files that have been read. 
  */  
  double * Extinction = (double *)calloc( No_of_SPECIES,sizeof(double));
  double * Colonization = (double *)calloc(No_of_SPECIES,sizeof(double));
  double * NLL_Value = (double *)calloc( No_of_SPECIES, sizeof(double));
  
  double Total_NLL_Value = 0.0;	
  for( i=0; i<No_of_SPECIES; i++ ){
    /* Initial colonization and extinction values are used to seed 
       the heuristic search. This values are given as input arguments
    */
    Colonization[i] = Colonization_Rate;      
    Extinction[i]   = Extinction_Rate;
    
    mle_Col_Ext_Uneven_Matrix_R_SHLIB( Presence[i], No_of_SITES[i],
				       Time_Vector[i], No_of_TIMES[i],
				       Sp_Time[i], No_of_Sp_Times[i], 
				       &Colonization[i], C_Range, 
				       &Extinction[i], E_Range,
				       No_of_PARAMETERS, No_of_PARAMETERS_MAX,
				       Index, Discretization, 
				       Tolerance, No_of_ITERATIONS, 
				       Verbose, Minimization, 
				       &NLL_Value[i] ); 
    Total_NLL_Value += NLL_Value[i];	           
    if( (*Verbose) == 1) {
      Rprintf(" File %d (%s): ", i, Species_Tag[i]);
      Rprintf(" NLL (Colonization = %g, Extinction = %g) = %g\n",  
	     Colonization[i], Extinction[i], NLL_Value[i] ); 
    }
  }
  /*     E N D : --------------------------------------------------------
   */
  
  /* B E G I N : --------------------------------------------------------
     3. Calculation of a between-species distance matrix on the 
     colonization-extinction parameter space   
  */
  double ** Distance_Matrix  = (double **)calloc( No_of_SPECIES, 
						  sizeof(double *) );
  for(i=1; i<No_of_SPECIES; i++) 
    Distance_Matrix[i] = (double *)calloc( i, sizeof(double) );
  
  EUCLIDIAN_DISTANCE_MATRIX( Colonization, Extinction, No_of_SPECIES,
			     Distance_Matrix );
  /*     E N D : ------------------------------------------------------- 
   */  
  
  /* B E G I N : -------------------------------------------------------
     4. The distance matrix is used to feed a typical UPGMA clustering 
     algorithm and the binary grouping is then used to generate models 
     of different number of parameters and, finally, the functions 
     conducts standard AIC-based model selection.
  */
  MODEL_SELECTION_UPGMA_R_SHLIB( Colonization_Rate, C_Range, 
				 Extinction_Rate, E_Range,  
				 Data, Distance_Matrix, No_of_SPECIES,
				 No_of_PARAMETERS, No_of_PARAMETERS_MAX,
				 Index, Discretization, 
				 Tolerance, No_of_ITERATIONS, 
				 Verbose, Minimization,
				 Results ); 
  /*    E N D : --------------------------------------------------------
   */
  
  /* B E G I N : -------------------------------------------------------
     Freeing up allocated memmory!!! 
  */
  for( i=0; i<No_of_SPECIES; i++ ){
    free(No_of_Sp_Times[i]);
    for (j=0; j<No_of_SITES[i]; j++) free( Sp_Time[i][j] );
    free( Sp_Time[i] );
  }
  free(Sp_Time); free(No_of_Sp_Times);
  
  free(Colonization); 
  free(Extinction); 
  free(NLL_Value);
  
  for(i=1; i<No_of_SPECIES; i++) free(Distance_Matrix[i]);
  free(Distance_Matrix);
  
  for( i=0; i<No_of_SPECIES; i++ ) SP_Matrix_Data_Free( Data[i] ); 
  free (Data); 
  free( Dummy );
  /*   END: ------------------------------------------------------------ 
   */
}
