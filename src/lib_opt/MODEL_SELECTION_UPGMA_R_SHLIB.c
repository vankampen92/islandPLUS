#include <MODEL.h>

void MODEL_SELECTION_UPGMA_R_SHLIB( double Colonization_Rate, double * C_Range, 
				    double Extinction_Rate, double * E_Range, 
				    SP_Matrix_Data ** Data,			      
				    double ** Distance_Matrix,
				    int No_of_SPECIES, 
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
  /* This function calculates the MLE (Extinction, 
     Colonization) pair for each group of species and subgroupings
     Associated AIC statistics are also calculated. The partition is 
     defined by the tree that results from the UPGMA
     algortithm from the distance matrix 'Distance Matrix[][]'.
     
     Input Parameters:
      . Colonization_Rate, initial seed for the searches.       
      . Extinction_Rate, initial seed for the searches.
      . No_of_SPECIES, nombre of different basic objects to classify. 
        This is also, the total number of different partitions in which 
        a set of different objects can be grouped if subdivisions are 
        made in a binary way. In this case, we know that the total 
	number of possible partitions matches the total number of 
	initial elements to be partitioned (No_of_SPECIES).   
      . Distance_Matrix[][]
      . Data, an array of data structures of lenght No_of_SPECIES
        that simplifies communication between procedures. 

     Output: 
      . A file, in LaTex format, is produced containing a table with  
        NLL, AIC, AIC_c, AIC_d differences and AIC_w model weights 
	in columns where each row corresponds to a different model, 
	this is, a different partition, with a different number of 
	parameters. Each partition defines a different model with a 
	different number of parameters. Akaike Information Criteria 
	are used for model comparison.
      . Results[][]: A matrix containing in rows (with no columns 
        names) 6 columns: No of Parameters NLL AIC AIC_c AIC_d AIC_w
  
     The function is organized in three steps:

     1. The distance matrix is used to feed a typical UPGMA clustering
     algorithm. Importantly, the output of this procedure characterizes
     a partition (or binary tree) by means of 3 arrays (PARTITION[][][], 
     G[], and K[][], see below).

     2.  This procedure calculates the MLE (Extinction, Colonization) 
     pair for each group of species for each and every partition as 
     defined by PARTITION above. Each partition is then characterized
     by its NLL, AIC, AIC_c, AIC_d differences and AIC_w weights. 
     In particular, an index i loops over partitions (from i=0 to No_of_SPECIES-1). 
     Within this big loop:

        a. A total Presence matrix is generated lumping together individual 
	   Presence matrices. Partition_Presence[i][j] is the matrix of the j-th group 
	   (which may include several species) belonging to the i-th partition. 

	b. Partition_Presence[i][j] is passed over the mle function:
	   mle_Col_Ext_Uneven_Matrix( Partition_Presence[i][j], R[i][j],
	                              General_Time_Vector, 
	                 	      General_No_of_SAMPLING_TIMES,
				      Sp_Time_Vector[i][j], No_of_Sp_Times[i][j], 
				      &Colonization[i][j], &Extinction[i][j], 
				      &NLL_Value[i][j] );

       As a consequence, as the output of this main loop, the MLE 
       (colonization, extinction) pair is calculated for each group for every 
       partition and its corresponding NLL total value. 
       
      3. Finally, NLL, AIC, AIC_c, AIC_d differences and AIC_w model weights are 
      organized in a table (each row corresponds to a different model, this is, to a 
      different partition, with a different number of parameters), which 
      is saved in a file (latex format). 
  */ 
  int i,j,k,l,s,  m,n;
  /* B E G I N : -----------------------------------------------------
     1. The distance matrix is used to feed a usual UPGMA clustering
     algorithm. Importantly, the output of this procedure is a 
     a full binary partition or tree, which is characterized by:

     . PARTITION[][][], where, for instance, PARTITION[i][j][k] is 
       the label of the k-th species belonging to the j-th subgroup of 
       the i-th partition. 
     . G[], where, G[i] is the total number of groups of the i-th 
       partition. For instance, if G[i]=2, it means the i-th partition 
       is made up by two subgroups. 
     . K[][], where K[i][j] is the number of single elements or the 
       cardinal of the j-th subgroup belonging the i-th partititon. 
  */
  int * G   = (int *)calloc( No_of_SPECIES, sizeof(int) );
  int ** K          = (int **)calloc( No_of_SPECIES, sizeof(int *) );
  int *** PARTITION = (int ***)calloc( No_of_SPECIES, sizeof(int **) );
  for (i=0; i<No_of_SPECIES; i++ ) {
    G[i] = i+1;
    K[i] = (int *)calloc( G[i], sizeof(int) );
    PARTITION[i] = (int **)calloc( G[i], sizeof(int *) );
    for (j=0; j<G[i]; j++ )
      PARTITION[i][j] = (int *)calloc( No_of_SPECIES, sizeof(int) );
  }
  
  UPGMA_CLUSTERING_PARTITION( Distance_Matrix, No_of_SPECIES,
			      PARTITION, G, K );
  for(i=0; i<No_of_SPECIES; i++) {
    if( (*Verbose) == 1 ) Rprintf(" %dth Partition is made up by %d subsets: ", 
				 i, G[i]);
    for(j=0; j<G[i]; j++) {
      if( (*Verbose) == 1 ) Rprintf("{ ");
      for(k=0; k<K[i][j]; k++) { 
	n = PARTITION[i][j][k]; 
	if( (*Verbose) == 1 ) Rprintf("%s ", Data[n]->Name);
      }
      if( (*Verbose) == 1 ) Rprintf("} ");
    }
    if( (*Verbose) == 1 ) Rprintf("\n");
  }
  /*     E N D : -----------------------------------------------------
   */

  /* B E G I N : -----------------------------------------------------
      2. Reserving memory and preparing it for the main partition loop (see below) 
  */
  /* B E G I N: Reserving memory ---------------------- */		     		 
  double ** NLL          = (double **)calloc( No_of_SPECIES, sizeof(double *) );
  double ** COL          = (double **)calloc( No_of_SPECIES, sizeof(double *) );
  double ** EXT          = (double **)calloc( No_of_SPECIES, sizeof(double *) );
  double * MODEL_NLL     = (double *)calloc( No_of_SPECIES, sizeof(double) );
  double * MODEL_AIC     = (double *)calloc( No_of_SPECIES, sizeof(double) );
  double * MODEL_AIC_c   = (double *)calloc( No_of_SPECIES, sizeof(double) );
  for (i=0; i<No_of_SPECIES; i++ ) {
    NLL[i] = (double *)calloc( G[i], sizeof(double) );
    COL[i] = (double *)calloc( G[i], sizeof(double) );
    EXT[i] = (double *)calloc( G[i], sizeof(double) );
  }
  int General_No_of_SAMPLING_TIMES = Data[0]->No_of_TIMES; 
  /* For instance, 17 different sampling times for zooplancton data. 
     All matrices should have the same number of sampling times (columns) 
     before filtering out missing values. See assert() below just in 
     case. 
  */
  for(i=0; i<No_of_SPECIES; i++) 
    if( General_No_of_SAMPLING_TIMES != Data[i]->No_of_TIMES ) error(0,0,"Number of Sampling Times do not match: program aborted");
  /* Calculating total number of rows of the presence matrix corresponding 
     to the j-th group of the i-th partition: R[i][j] */
  int ** R          = (int **)calloc( No_of_SPECIES, sizeof(int *) );
  for (i=0; i<No_of_SPECIES; i++ ) R[i] = (int *)calloc( G[i], sizeof(int) );
  int N_SUM;
  for( i=0; i<No_of_SPECIES; i++ ) {
    for( j=0; j<G[i]; j++ ) {
      N_SUM = 0;
      for( k=0; k<K[i][j]; k++ ) {
	n = PARTITION[i][j][k];
	N_SUM += Data[n]->No_of_SITES;
      }
      R[i][j] = N_SUM;
    }
  }
  
  double **** Partition_Presence = (double ****)calloc( No_of_SPECIES, 
							sizeof(double ***) );
  for( i=0; i<No_of_SPECIES; i++ ){
    Partition_Presence[i] = (double ***)calloc(G[i], sizeof(double **) );
    for( j=0; j<G[i]; j++ ) {
      Partition_Presence[i][j] = (double **)calloc( R[i][j], sizeof(double *));
      for( k=0; k<R[i][j]; k++ ){
	Partition_Presence[i][j][k] = (double *)calloc(General_No_of_SAMPLING_TIMES, 
						       sizeof(double ));	  
      } 
    }
  }
  /*     E N D: Reserving memory ---------------------- */		     

  /* B E G I N : --------------------------------------------------------  
     Partition (i) loop (the core of this function) 
  */
    for(i=0; i<No_of_SPECIES; i++) { /* i-th Partitition: 
      Since partition are made by a upwma algorithm, 
      the number of different partitions matches the number of species!!! */
      
      /* B E G I N: Reserving memory ---------------------- */		     
      double * General_Time_Vector = (double *)calloc( General_No_of_SAMPLING_TIMES, 
						       sizeof(double) ); 
      for (j=0; j<General_No_of_SAMPLING_TIMES; j++ ) {
	General_Time_Vector[j] = Data[i]->Time_Vector[j];
      }
      //General Time Vector should be the same vector across all species.
      //Warning: this is not checked here! 
      //If this is not true, the whole procedure is not reliable. 
      //Notice that the i-th index here stands for the i-th partition, which 
      //will be made up by a number of groups of different species. All these 
      //species should share the same Time Vector. If this is not true, 
      //some previous rearrangement of the input data file is compulsory. 
      
      double *** Sp_Time_Vector = (double ***)calloc(G[i], sizeof(double **));
      for(j=0; j<G[i]; j++) {
	Sp_Time_Vector[j] = (double **)calloc(R[i][j], sizeof(double *) );
	for (k=0; k<R[i][j]; k++) 
	  Sp_Time_Vector[j][k] = (double *)calloc(General_No_of_SAMPLING_TIMES, 
						  sizeof(double ));
      }
      
      int ** No_of_Sp_Times = (int **)calloc(G[i], sizeof(int *));
      for( j=0; j<G[i]; j++ )
	No_of_Sp_Times[j] = (int *)calloc(R[i][j], sizeof(int) );
      /*     E N D: Reserving memory ---------------------- */		     
      
      /* B E G I N : ----------------------------------------------------- 
	 Building Partition_Matrices associated to each group for every 
	 partition. In particular, here, associated to the i-th partition. 
      */
      if( (*Verbose) == 1 ) Rprintf(" %dth Partition is made up by %d subsets: ", 
				   i, G[i]);
      for(j=0; j<G[i]; j++) {
	if( (*Verbose) == 1 ) Rprintf("{ ");
	for(k=0; k<K[i][j]; k++) 
	  if( (*Verbose) == 1 ) Rprintf("%d ", PARTITION[i][j][k]);
	if( (*Verbose) == 1 ) Rprintf("} ");
	m = 0;
	for(k=0; k<K[i][j]; k++) {
	  n = PARTITION[i][j][k];  
	  for(l=0; l<Data[n]->No_of_SITES; l++) {
	    for(s=0; s<Data[n]->No_of_TIMES; s++) {  
	      Partition_Presence[i][j][m][s] = Data[n]->Presence[l][s];
	    }
	    m++;
	  }
	}
	if ( R[i][j] != m ) error(0,0,"Program aborted");
      }
      if( (*Verbose) == 1) Rprintf("\n");
      /*     E N D : -------------------------------------------------------- 
       */
      
      /* B E G I N : Dealing with missing data. Calculating true Time Vectors */
      int No_of_TRANSITIONS;
      int Total_No_of_TRANSITIONS;
      int N00, N01, N10, N11;
      Total_No_of_TRANSITIONS = 0;
      for( j=0; j<G[i]; j++ ) {
	IO_Filtering_Out_Matrix( Partition_Presence[i][j], &R[i][j], 
				 General_Time_Vector, General_No_of_SAMPLING_TIMES,  
				 Sp_Time_Vector[j],  No_of_Sp_Times[j], 
				 0.1 );
	No_of_TRANSITIONS = 0;
	for( k=0; k<R[i][j]; k++ ) {
	  
	  if(No_of_Sp_Times[j][k] <= 1) error(0,0,"Program aborted");

	  if( (*Verbose) == 1 ) Rprintf(" Times:\t"); 
	  for(n=0; n<No_of_Sp_Times[j][k]; n++) 
	    if( (*Verbose) == 1 ) Rprintf("%g ", Sp_Time_Vector[j][k][n]);
	  if( (*Verbose) == 1 ) Rprintf("\n");
	  if( (*Verbose) == 1 ) Rprintf(" 0 / 1:\t"); 
	  for(n=0; n<No_of_Sp_Times[j][k]; n++) 
	    if( (*Verbose) == 1 ) Rprintf("%g ", Partition_Presence[i][j][k][n]);
	  if( (*Verbose) == 1 ) Rprintf("\n");

	  No_of_TRANSITIONS += (No_of_Sp_Times[j][k]-1);
	}
	  Total_No_of_TRANSITIONS += No_of_TRANSITIONS; 
	  
	  Counting_Type_of_Transitions( Partition_Presence[i][j], R[i][j], 
	    Sp_Time_Vector[j], No_of_Sp_Times[j], 
	    &N00, &N01, &N10, &N11 );
	  
	  if( N00 + N01 + N10 + N11 != No_of_TRANSITIONS ) error(0,0,"Number of Transititions do not sum up");
	  
	  if( (*Verbose) == 1 ) 
	    Rprintf(" Total No of TRANSITIONS (Partition: %d-th: Group %d-th) = %d\n", 
	    i, j, No_of_TRANSITIONS );
	  if( (*Verbose) == 1 ) 
	    Rprintf(" No of observed extinctions       ( 1 ---> 0 ): %d\n", N01);
	  if( (*Verbose) == 1 ) 
	    Rprintf(" No of observed colonizations     ( 0 ---> 1 ): %d\n", N10);
	  if( (*Verbose) == 1 ) 
	    Rprintf(" No of observed non-colonizations ( 0 ---> 0 ): %d\n", N00);
	  if( (*Verbose) == 1 ) 
	    Rprintf(" No of observed permanences       ( 1 ---> 1 ): %d\n", N11);
	  if( (*Verbose) == 1 ) 
	    Rprintf("\n\n");
	  if( (*Verbose) == 1 ) 
	    getchar();
	}
      /*     E N D : ---------------------------------------------------------*/
      
      /* B E G I N : Calculating mle by using General_No_of_SAMPLING_TIMES columns */ 
      MODEL_NLL[i] = 0.0;	
      for( j=0; j<G[i]; j++ ){
	/* B E G I N : ------------------------------------------ 
	   Initial colonization and extinction values are used to
	 seed the heuristic search 
	*/
	COL[i][j] = Colonization_Rate;      
	EXT[i][j] = Extinction_Rate;
	/*     E N D : ------------------------------------------ 
	 */ 
	mle_Col_Ext_Uneven_Matrix_R_SHLIB( Partition_Presence[i][j], R[i][j],
					   General_Time_Vector, 
					   General_No_of_SAMPLING_TIMES,  
					   Sp_Time_Vector[j], No_of_Sp_Times[j], 
					   &COL[i][j], C_Range,  
					   &EXT[i][j], E_Range,
					   No_of_PARAMETERS, No_of_PARAMETERS_MAX,
					   Index, Discretization, 
					   Tolerance, No_of_ITERATIONS, 
					   Verbose, Minimization, 
					   &NLL[i][j] );
	MODEL_NLL[i] += NLL[i][j];	           
	if( (*Verbose) == 1) Rprintf("Partition: %d-th: Group %d-th:", i,j );
	if( (*Verbose) == 1) Rprintf(" NLL (Colonization = %g, Extinction = %g) = %g\n",  
				    COL[i][j], EXT[i][j], NLL[i][j] ); 
      }
      /*     E N D : End Index Calculation                                     */
      
      if( (*Verbose) == 1 ) Rprintf(" %dth Partition (%d subsets): ", i, G[i]);
      for(j=0; j<G[i]; j++) { 
	if( (*Verbose) == 1 ) Rprintf("{ ");
	for(k=0; k<K[i][j]; k++) if( (*Verbose) == 1 ) Rprintf("%d ", PARTITION[i][j][k]);
	if( (*Verbose) == 1 ) Rprintf("} ");
      }
      if( (*Verbose) == 1 ) Rprintf("\n");
      if( (*Verbose) == 1 ) Rprintf(" The total negative loglikelihood of partition %dth:\n", i);
      if( (*Verbose) == 1 ) Rprintf(" NLogL(M | %d-th Partition) = %g\n", 
				   i, MODEL_NLL[i]);
      if( (*Verbose) == 1 ) Rprintf(" The total number of transtions, i.e., the total number of factors in the multiplicative likelihood is: %d\n", 
				   Total_No_of_TRANSITIONS);
      
      int No_of_Estimated_Parameters = G[i] * 2;
      double Ka = (double)No_of_Estimated_Parameters;
      MODEL_AIC[i] = 2.0 * MODEL_NLL[i] + 2.0 * Ka;
      MODEL_AIC_c[i] = MODEL_AIC[i] +  2.0 * Ka * (Ka + 1) / ((double)Total_No_of_TRANSITIONS - Ka - 1.0);
      if( (*Verbose) == 1 ) Rprintf(" Partition %d-th: Number of estimated parameters: %d\n AIC = %g\tAIC (corrected) = %g\n", 
				   i, No_of_Estimated_Parameters, 
				   MODEL_AIC[i], MODEL_AIC_c[i]);
      //getchar();

      for( j=0; j<G[i]; j++ ){
	for (k=0; k<R[i][j]; k++) free(Sp_Time_Vector[j][k]);
	free(Sp_Time_Vector[j]);
      }
      free (Sp_Time_Vector);

      for( j=0; j<G[i]; j++ ) free (No_of_Sp_Times[j]);
      free (No_of_Sp_Times);

      free( General_Time_Vector );
    }
    /*     E N D : Partition (i) loop * * * * * * * * * * * * * * * * * * * * * */
    
    char ** Name = (char **)calloc( No_of_SPECIES, sizeof(char *) );
    for (i=0; i<No_of_SPECIES; i++) {
      Name[i] = (char *)calloc( 5, sizeof(char) );
      memcpy( Name[i], Data[i]->Name, strlen(Data[i]->Name) );
    }
    Model_Selection_AIC_Latex_Table( Name, 
				     PARTITION, G, K, No_of_SPECIES, 
				     MODEL_NLL, COL, EXT, 
				     MODEL_AIC, MODEL_AIC_c ); 
    
    int * No_of_MODEL_PARAMETERS = (int *)calloc( No_of_SPECIES, sizeof(int) );
    for(i=0; i<No_of_SPECIES; i++) { 
      No_of_MODEL_PARAMETERS[i] = 2*G[i]; 
    }
    AIC_Summary_RESULTS( No_of_MODEL_PARAMETERS, No_of_SPECIES,  
			 MODEL_NLL, MODEL_AIC, MODEL_AIC_c, 
			 Results );  

    /* B E G I N: Freeing up memory ---------------------- */
    free(No_of_MODEL_PARAMETERS); 
    for (i=0; i<No_of_SPECIES; i++) free(Name[i]); 
    free(Name);
    
    for( i=0; i<No_of_SPECIES; i++ ){
      for( j=0; j<G[i]; j++ ) {
	for( k=0; k<R[i][j]; k++ ) free( Partition_Presence[i][j][k] );
	free(Partition_Presence[i][j]);
      }
      free( Partition_Presence[i] );
    }
    free( Partition_Presence ); 

    for(i = 0; i<No_of_SPECIES; i++) {
      free(NLL[i]);
      free(COL[i]);
      free(EXT[i]);
      free(R[i]);
    }
    free(NLL); free(COL); free(EXT); free(R);
    free(MODEL_NLL);  
    free(MODEL_AIC);  
    free(MODEL_AIC_c);			       
    /*   E N D: ---------------------------------------- */
}
