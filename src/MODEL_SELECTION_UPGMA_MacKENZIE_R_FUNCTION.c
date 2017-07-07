#include <MODEL.h>

void MODEL_SELECTION_UPGMA_MacKENZIE_DRIVER ( int No_of_SPECIES, char ** Species_Tag, 
					      double *** Presence, int * No_of_SITES, 
					      double ** Time_Vector, int * No_of_TIMES, 
					      double MISSING_VALUE_FLAG, 
					      double Colonization_Rate,
					      double * C_Range, 
					      double Extinction_Rate,
					      double * E_Range,
					      double Detectability,
					      double * D_Range, 
					      double Phi_Time_0,
					      double * P_Range,
					      int * No_of_PARAMETERS,
					      int * No_of_PARAMETERS_MAX,
					      int * Index, 
					      int * Discretization,
					      double * Tolerance, 
					      int * No_of_ITERATIONS,
					      int * Verbose, 
					      int * Minimization,
					      double ** Results );

void MODEL_SELECTION_UPGMA_MacKENZIE_R_FUNCTION ( int * pNo_of_SPECIES,
						  char ** Species_Tag, 
						  double * R_Presence,
						  int * No_of_SITES, 
						  double * R_Time_Vector,
						  int * No_of_TIMES, 
						  double * MISSING_VALUE_FLAG, 
						  double * Colonization_Rate,
						  double * C_Range,
						  double * Extinction_Rate,
						  double * E_Range,
						  double * Detectability,
						  double * D_Range, 
						  double * Phi_Time_0,
						  double * P_Range,
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
  /* Input parameters:
     . No_of_SPECIES, this is the number of elemental/single groups considered 
     . Species_Tag, this the identifier of each and every of these single groups

     (the rest of input parameters are self-explanatory)
  */
  
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
    No_of_COLUMNS = MAX( No_of_COLUMNS, No_of_TIMES[i] );  
  
  Rprintf( "Number of Columns: %d\n", No_of_COLUMNS );  //getchar();

  n=0; 
  for(i=0; i<No_of_SPECIES; i++) 
    for(j=0; j<No_of_SITES[i]; j++)
      for(k=0; k<No_of_TIMES[i]; k++)
	Presence[i][j][k] = R_Presence[n++];

  MODEL_SELECTION_UPGMA_MacKENZIE_DRIVER (No_of_SPECIES, Species_Tag, 
					  Presence, No_of_SITES, 
					  Time_Vector, No_of_TIMES, 
					  (* MISSING_VALUE_FLAG), 
					  (* Colonization_Rate), C_Range,  
					  (* Extinction_Rate), E_Range,
					  (* Detectability), D_Range, 
					  (* Phi_Time_0), P_Range,
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

void MODEL_SELECTION_UPGMA_MacKENZIE_DRIVER (int No_of_SPECIES, char ** Species_Tag, 
					     double *** Presence, int * No_of_SITES, 
					     double ** Time_Vector, int * No_of_TIMES, 
					     double MISSING_VALUE_FLAG, 
					     double Colonization_Rate,
					     double * C_Range, 
					     double Extinction_Rate,
					     double * E_Range,
					     double Detect_Probability,
					     double * D_Range, 
					     double Phi_Time_0,
					     double * P_Range,
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
  int ** Sp_No_of_Times = (int **)calloc(No_of_SPECIES, sizeof(int *));
  for( i=0; i<No_of_SPECIES; i++ ) {
    Sp_No_of_Times[i] = (int *)calloc(No_of_SITES[i], sizeof(int) );
    for( j=0; j<No_of_SITES[i]; j++) {
      Sp_No_of_Times[i][j] = No_of_TIMES[i];
    }
  }

  int ** Transects = (int **)calloc(No_of_SPECIES, sizeof(int *));
  for( i=0; i<No_of_SPECIES; i++ ) {
    Transects[i] = (int *)calloc(No_of_TIMES[i], sizeof(int) );
  }
  
  int *** Sp_Transect = (int ***)calloc(No_of_SPECIES, sizeof(int **));
  for( i=0; i<No_of_SPECIES; i++ ){
    Sp_Transect[i] = (int **)calloc(No_of_SITES[i], sizeof(int *) );
    for (j=0; j<No_of_SITES[i]; j++) { 
      Sp_Transect[i][j] = (int *)calloc(No_of_TIMES[i], sizeof(int));
    }
  }

  int ** Sp_No_of_Transects = (int **)calloc(No_of_SPECIES, sizeof(int *));
  for( i=0; i<No_of_SPECIES; i++ ) {
    Sp_No_of_Transects[i] = (int *)calloc(No_of_SITES[i], sizeof(int) );
  }

  /* B E G I N : -------------------------------------------------------
     1: Calculating true time vectors and true transect vectors per row 
        of each of the data matrices for each group. 
   
	Before filtering out missing values, full data matrices and sampling 
	times are allocated in an array of length No_of_SPECIES of type 
	SP_Matrix_Data (see definition in MODEL_SP_Matrix_Data_STRUCT.h). This  
	generic data structure allows to store all these data information in 
	a convenient manner. The reason for that is that, after filetering out 
	missing values, Presence matrices are changed on place. As all these
	operations are done, after this section, fully initialization of an 
	array of data structures named 'Data' of size No_of_SPECIES, this is, 
	the number of single groups, is achieved. 
  */
  SP_Matrix_Data ** Data = (SP_Matrix_Data **)calloc( No_of_SPECIES, 
						      sizeof(SP_Matrix_Data *) );
  int No_of_COLUMNS = No_of_TIMES[0]; 
  for(i=0; i<No_of_SPECIES; i++)
    No_of_COLUMNS = MAX( No_of_COLUMNS, No_of_TIMES[i] );

  int * Total_No_of_Transects  = (int *)calloc( No_of_COLUMNS, sizeof(int) );
  for(i=0; i<No_of_SPECIES; i++) Total_No_of_Transects[i] = No_of_TIMES[i]; 
  
  double * Dummy = (double *)calloc( No_of_COLUMNS, sizeof(double) );
  int    Dummy_No_of_TIMES = No_of_COLUMNS; 
  for(i=0; i<No_of_SPECIES; i++) {
    for(j=0; j<No_of_TIMES[i]; j++) Dummy[j] = Time_Vector[i][j];
      
    Counting_Replicates_per_Time(Time_Vector[i], &No_of_TIMES[i],
				 Dummy, Transects[i], 
				 &Dummy_No_of_TIMES);
  }
  
  for( i=0; i<No_of_SPECIES; i++ ) {
    Data[i] = SP_Matrix_Data_Alloc( No_of_SITES[i], No_of_COLUMNS,
				    No_of_COLUMNS ) ;
    SP_Matrix_Data_Setup( No_of_SITES[i], No_of_COLUMNS, Total_No_of_Transects[i], 
			  Data[i], Presence[i], 
			  Time_Vector[i], Sp_Time[i], Sp_No_of_Times[i], 
			  Transects[i], Species_Tag[i] );
  }
  
  IO_Filtering_Out_Missing_Values ( No_of_SPECIES,  
				    Presence, No_of_SITES, 
				    Time_Vector, No_of_TIMES, 
				    Sp_Time, Sp_No_of_Times, 
				    MISSING_VALUE_FLAG );

  for(i=0; i<No_of_SPECIES; i++) {
    for(j=0; j<No_of_SITES[i]; j++) {
      Sp_No_of_Transects[i][j] = Sp_No_of_Times[i][j]; 
      Counting_Replicates_per_Time(Sp_Time[i][j], &Sp_No_of_Transects[i][j],
				   Sp_Time[i][j], Sp_Transect[i][j], 
				   &Sp_No_of_Times[i][j] );
    }
  }

  for( i=0; i<No_of_SPECIES; i++ ) {
    
    Total_No_of_Transects[i] = No_of_TIMES[i]; 
    Counting_Replicates_per_Time(Time_Vector[i], &Total_No_of_Transects[i],
				 Time_Vector[i], Transects[i], 
				 &No_of_TIMES[i]);
    
    SP_Matrix_Data_Uneven_Setup( Data[i], Presence[i], 
				 No_of_SITES[i], Total_No_of_Transects[i], 
				 Time_Vector[i], Transects[i], No_of_TIMES[i],
				 Sp_Time[i], Sp_No_of_Times[i],
				 Sp_Transect[i], Sp_No_of_Transects[i] );

    if( (* Verbose) == 1 ) SP_Matrix_Data_Write( Data[i] );

    Rprintf(" Data structure for the %d-th elemental group (%s) is done\n\n",
	    i, Species_Tag[i] );
  }
  /*       E N D :------------------------------------------------------
   */

  /* B E G I N : -------------------------------------------------------
     2: Calculation of mle of extinction, colonization, detectability and 
        phi_0 model parameters: four parameters per each of the groups that 
	have been considered. 
  */  
  double * Extinction = (double *)calloc( No_of_SPECIES,sizeof(double));
  double * Colonization = (double *)calloc(No_of_SPECIES,sizeof(double));
  double * Detectability = (double *)calloc( No_of_SPECIES,sizeof(double));
  double * Phi_T_0 = (double *)calloc(No_of_SPECIES,sizeof(double));
  double * NLL_Value = (double *)calloc( No_of_SPECIES, sizeof(double));
  
  double Total_NLL_Value = 0.0;	
  for( i=0; i<No_of_SPECIES; i++ ){
    /* Initial colonization and extinction values are used to seed 
       the heuristic search. This values are given as input arguments
    */
    Colonization[i]  = Colonization_Rate;      
    Extinction[i]    = Extinction_Rate;
    Detectability[i] = Detect_Probability;
    Phi_T_0[i]    = Phi_Time_0; 
    
    mle_MacKenzie_Uneven_Matrix_R_SHLIB( Presence[i],
					 No_of_SITES[i], Total_No_of_Transects[i], 
					 Time_Vector[i], Transects[i], No_of_TIMES[i],
					 Sp_Time[i], Sp_No_of_Times[i],
					 Sp_Transect[i], Sp_No_of_Transects[i],  
					 &Colonization[i], C_Range, 
					 &Extinction[i], E_Range,
					 &Detectability[i], D_Range, 
					 &Phi_T_0[i], P_Range, 
					 No_of_PARAMETERS, No_of_PARAMETERS_MAX,
					 Index, Discretization, 
					 Tolerance, No_of_ITERATIONS, 
					 Verbose, Minimization, 
					 &NLL_Value[i] );
    
    Total_NLL_Value += NLL_Value[i];	           
    // if( (*Verbose) == 1) {
      Rprintf(" Group %d (%s): ", i, Species_Tag[i]);
      Rprintf(" NLL (Col = %g, Ext = %g, Dtc = %g, P_0 = %g) = %g\n",  
	      Colonization[i], Extinction[i], Detectability[i], Phi_T_0[i],
	      NLL_Value[i] ); 
    // }
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
  
  EUCLIDIAN_DISTANCE_MATRIX_MacKENZIE( Colonization, Extinction,
				       Detectability, Phi_T_0,
				       No_of_SPECIES, Distance_Matrix );
  /*     E N D : ------------------------------------------------------- 
   */  
  
  /* B E G I N : -------------------------------------------------------
     4. The distance matrix is used to feed a typical UPGMA clustering 
     algorithm and the binary grouping is then used to generate models 
     of different number of parameters and, finally, this function below 
     conducts standard AIC-based model selection.
  */
  MODEL_SELECTION_UPGMA_MacKENZIE_R_SHLIB( Colonization_Rate, C_Range, 
					   Extinction_Rate, E_Range,
					   Detect_Probability, D_Range, 
					   Phi_Time_0, P_Range, 
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
    free(Sp_No_of_Times[i]);
    free(Transects[i]);
    free(Sp_No_of_Transects[i]); 
    for (j=0; j<No_of_SITES[i]; j++) {
      free( Sp_Time[i][j] );
      free( Sp_Transect[i][j] );
    }
    free( Sp_Time[i] );
    free( Sp_Transect[i] );
  }
  free(Sp_Time); free(Sp_No_of_Times);
  free(Transects); free(Sp_No_of_Transects); free(Sp_Transect); 
  
  free(Colonization); 
  free(Extinction);
  free(Detectability);
  free(Phi_T_0);
  free(NLL_Value);
  
  for(i=1; i<No_of_SPECIES; i++) free(Distance_Matrix[i]);
  free(Distance_Matrix);
  
  for( i=0; i<No_of_SPECIES; i++ ) SP_Matrix_Data_Free( Data[i] ); 
  free (Data);

  free (Dummy);
  free (Total_No_of_Transects); 
  /*   END: ------------------------------------------------------------ 
   */
}
