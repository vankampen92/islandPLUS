#include <MODEL.h>

void MODEL_SELECTION_UPGMA_MacKENZIE_R_SHLIB(double Colonization_Rate, double * C_Range, 
					     double Extinction_Rate, double * E_Range,
					     double Detectability, double * D_Range, 
					     double Phi_Time_0, double * P_Range,
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
      . Detectability, initial seed for the searches.
      . Phi_Time_0, initial seed for the searches.
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
	   (which may include several elemental grops) belonging to the i-th partition. 

	b. Partition_Presence[i][j] is passed over the mle function:
	   mle_MacKenzie_Uneven_Matrix( Partition_Presence[i][j], R[i][j],
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
     different partition, with a different number of parameters), which is saved in 
     a file (latex format). 
  */ 
  int i,j,k,l,s,  m,n;
  int Total_No_of_TRANSITIONS;
  
  /* B E G I N : -----------------------------------------------------
     1. The distance matrix is used to feed a usual UPGMA clustering
     algorithm. Importantly, the output of this procedure is a 
     a full tree binary partition, which is characterized by:

     . PARTITION[][][], where, for instance, PARTITION[i][j][k] is 
       the label of the k-th species belonging to the j-th subgroup of 
       the i-th partition. 
     . G[], where, G[i] is the total number of groups of the i-th 
       partition. For instance, if G[i]=2, it means the i-th partition 
       is made up by two subgroups. 
     . K[][], where K[i][j] is the number of single elements or  
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
    if( (*Verbose) == 1 ) printf(" %dth Partition is made up by %d subsets: ", 
				 i, G[i]);
    for(j=0; j<G[i]; j++) {
      if( (*Verbose) == 1 ) printf("{ ");
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
  double ** DTC          = (double **)calloc( No_of_SPECIES, sizeof(double *) );
  double ** P_0          = (double **)calloc( No_of_SPECIES, sizeof(double *) );
  double * MODEL_NLL     = (double *)calloc( No_of_SPECIES, sizeof(double) );
  double * MODEL_AIC     = (double *)calloc( No_of_SPECIES, sizeof(double) );
  double * MODEL_AIC_c   = (double *)calloc( No_of_SPECIES, sizeof(double) );
  for (i=0; i<No_of_SPECIES; i++ ) {
    NLL[i] = (double *)calloc( G[i], sizeof(double) );
    COL[i] = (double *)calloc( G[i], sizeof(double) );
    EXT[i] = (double *)calloc( G[i], sizeof(double) );
    DTC[i] = (double *)calloc( G[i], sizeof(double) );
    P_0[i] = (double *)calloc( G[i], sizeof(double) );
  }
  int General_No_of_SAMPLING_TIMES = Data[0]->Total_No_of_TRANSECTS; 
  for(i=0; i<No_of_SPECIES; i++) 
    General_No_of_SAMPLING_TIMES = MAX(General_No_of_SAMPLING_TIMES,
				       Data[i]->Total_No_of_TRANSECTS);
  
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
   *  Partition (i) loop (the core of this function) 
   */
  for(i=0; i<No_of_SPECIES; i++) { /* i-th Partitition: 
					Since partition are made by a upwma algorithm, 
					the number of different partitions matches the 
					number of species (initial single groups)!!! */
      
    /* B E G I N: Reserving memory ---------------------- */		     
    double * General_Time_Vector = (double *)calloc( General_No_of_SAMPLING_TIMES, 
						     sizeof(double) ); 
    
    double *** Sp_Time_Vector = (double ***)calloc(G[i], sizeof(double **));
    for(j=0; j<G[i]; j++) {
      Sp_Time_Vector[j] = (double **)calloc(R[i][j], sizeof(double *) );
      for (k=0; k<R[i][j]; k++) {
	Sp_Time_Vector[j][k] = (double *)calloc(General_No_of_SAMPLING_TIMES, 
						sizeof(double));
      }
    }
    
    int ** Sp_No_of_Times = (int **)calloc(G[i], sizeof(int *));
    for( j=0; j<G[i]; j++ )
      Sp_No_of_Times[j] = (int *)calloc(R[i][j], sizeof(int) );
    
    int * Transects = (int *)calloc(General_No_of_SAMPLING_TIMES, sizeof(int));
    
    int *** Sp_Transect = (int ***)calloc(G[i], sizeof(int **));
    for( j=0; j<G[i]; j++ ){
      Sp_Transect[j] = (int **)calloc(R[i][j], sizeof(int *) );
      for (k=0; k<R[i][j]; k++) { 
	Sp_Transect[j][k] = (int *)calloc(General_No_of_SAMPLING_TIMES,
					  sizeof(int));
      }
    }
    
    int ** Sp_No_of_Transects = (int **)calloc(G[i], sizeof(int *));
    for( j=0; j<G[i]; j++ ) 
      Sp_No_of_Transects[j] = (int *)calloc(R[i][j], sizeof(int) );
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
	  for(s=0; s<Data[n]->Sp_Total_No_Transects[l]; s++) {  
	    Partition_Presence[i][j][m][s] = Data[n]->Presence[l][s];
	  }
	  m++;
	}
      }
      if( R[i][j] != m ) error(0,0,"Program aborted");
    }
    if( (*Verbose) == 1) Rprintf("\n"); 
    /*     E N D : -------------------------------------------------------- 
     */

    /* B E G I N : -------------------------------------------------------------------
     *             Counting the number of factors in the likelihood function of partit-
     *             ion i-th (Total_No_of_TRANSITIONS). 
     *             Initializing True Time Vectors and associated Transects Vectors 
     */
    Total_No_of_TRANSITIONS = 0;
    for( j=0; j<G[i]; j++ ) {
      m = 0; 
      for(k=0; k<K[i][j]; k++) {
	n = PARTITION[i][j][k];
	for(l=0; l<Data[n]->No_of_SITES; l++) {
	  Sp_No_of_Transects[j][m] = Data[n]->Sp_Total_No_Transects[l];
	  Sp_No_of_Times[j][m]     = Data[n]->No_Sp_Time[l];
	  for(s=0; s<Data[n]->No_Sp_Time[l]; s++) {  
	    Sp_Transect[j][m][s]   = Data[n]->Sp_Transects[l][s];
	    Sp_Time_Vector[j][m][s] = Data[n]->Sp_Time[l][s];
	    if(s>1)
	      Total_No_of_TRANSITIONS += Data[n]->Sp_Transects[l][s];
	  }
	  m++;
	}
	
      }
      assert( R[i][j] == m );
      
      if( (*Verbose) == 1) Rprintf("\n");
      m = 0; 
      for( k=0; k<K[i][j]; k++ ) {
	n = PARTITION[i][j][k];  
	for(l=0; l<Data[n]->No_of_SITES; l++) {
	  if( (*Verbose) == 1 ) Rprintf(" Times:\t"); 
	  for(s=0; s<Data[n]->No_Sp_Time[l]; s++){
	    if( (*Verbose) == 1 ) Rprintf("%g ", Sp_Time_Vector[j][m][s]);
	  }
	  if( (*Verbose) == 1 ) Rprintf("\n");
	  if( (*Verbose) == 1 ) Rprintf(" (1/0):\t"); 
	  for(s=0; s<Data[n]->Sp_Total_No_Transects[l]; s++) {  
	    if( (*Verbose) == 1 ) Rprintf("%g ", Partition_Presence[i][j][m][s]);
	  }
	  if( (*Verbose) == 1 ) Rprintf("\n");
	  m++; 
	}
      }
      if( R[i][j] != m ) error(0,0,"Program aborted"); // assert( R[i][j] == m );
      // if( (*Verbose) == 1 )  getchar();
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
      DTC[i][j] = Detectability;      
      P_0[i][j] = Phi_Time_0;
      /*     E N D : ------------------------------------------ 
       */ 
      mle_MacKenzie_Uneven_Matrix_R_SHLIB( Partition_Presence[i][j], R[i][j],
					   General_No_of_SAMPLING_TIMES,
					   General_Time_Vector, Transects, 
					   General_No_of_SAMPLING_TIMES,
					   Sp_Time_Vector[j], Sp_No_of_Times[j],
					   Sp_Transect[j], Sp_No_of_Transects[j], 
					   &COL[i][j], C_Range,  
					   &EXT[i][j], E_Range,
					   &DTC[i][j], D_Range,
					   &P_0[i][j], P_Range, 
					   No_of_PARAMETERS, No_of_PARAMETERS_MAX,
					   Index, Discretization, 
					   Tolerance, No_of_ITERATIONS, 
					   Verbose, Minimization, 
					   &NLL[i][j] );
      MODEL_NLL[i] += NLL[i][j];	           
      /* if( (*Verbose) == 1) Rprintf("Partition: %d-th: Group %d-th:", i,j );              */
      /* if( (*Verbose) == 1) Rprintf(" NLL (C = %g, E = %g, D = %g, P = %g) = %g\n",   */
      /* 				   COL[i][j], EXT[i][j], DTC[i][j], P_0[i][j], NLL[i][j] );               */
      Rprintf("Partition: %d-th: Group %d-th:", i,j );
      Rprintf(" NLL (C = %g, E = %g, D = %g, P = %g) = %g\n",
	      COL[i][j], EXT[i][j], DTC[i][j], P_0[i][j], NLL[i][j] ); 
    }
    
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
    
    int No_of_Estimated_Parameters = G[i] * (* No_of_PARAMETERS);
    double Ka = (double)No_of_Estimated_Parameters;
    MODEL_AIC[i] = 2.0 * MODEL_NLL[i] + 2.0 * Ka;
    MODEL_AIC_c[i] = MODEL_AIC[i] +  2.0 * Ka * (Ka + 1) / ((double)Total_No_of_TRANSITIONS - Ka - 1.0);
    if( (*Verbose) == 1 ) Rprintf(" Partition %d-th: Number of estimated parameters: %d\n AIC = %g\tAIC (corrected) = %g\n", 
				 i, No_of_Estimated_Parameters, 
				 MODEL_AIC[i], MODEL_AIC_c[i]);
    //getchar();
    /* B E G I N : -------------------------------------------------------
       Freeing up allocated memmory (associated to the i-th partitiion!!!) 
    */
    for( j=0; j<G[i]; j++ ){
      for (k=0; k<R[i][j]; k++) free(Sp_Time_Vector[j][k]);
      free(Sp_Time_Vector[j]);
    }
    free (Sp_Time_Vector);
    
    for( j=0; j<G[i]; j++ ) free (Sp_No_of_Times[j]);
    free (Sp_No_of_Times);
    
    free( General_Time_Vector );  free (Transects);
    
    for( j=0; j<G[i]; j++ ){
      for (k=0; k<R[i][j]; k++)  free(Sp_Transect[j][k]); 
      free(Sp_Transect[j]); 
    }
    free(Sp_Transect); 
  
    for( j=0; j<G[i]; j++ ) free(Sp_No_of_Transects[j]);
    free( Sp_No_of_Transects);
    /*   END: ------------------------------------------------------------ 
     */
  }
  /*     E N D : -----------------------------------------------------------
   *     Partition (i) loop (the core of this function) 
   */
    
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
    No_of_MODEL_PARAMETERS[i] =  (* No_of_PARAMETERS) * G[i]; 
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
    free(DTC[i]);
    free(P_0[i]);
    free(R[i]);
  }
  free(NLL); free(COL); free(EXT); free(DTC); free(P_0); free(R);
  free(MODEL_NLL);  
  free(MODEL_AIC);  
  free(MODEL_AIC_c);
  
  for (i=0; i<No_of_SPECIES; i++ ) {
    for (j=0; j<G[i]; j++ ) free (PARTITION[i][j]); 
    free( PARTITION[i] );
  }
  free(G); free(K); free(PARTITION); 
  /*   E N D: ---------------------------------------- */
}
