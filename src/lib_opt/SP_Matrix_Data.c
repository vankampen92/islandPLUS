#include <MODEL.h>

SP_Matrix_Data * SP_Matrix_Data_Alloc( int No_of_SITES, 
				       int No_of_TIMES,
				       int Total_No_of_TRANSECTS )
{
  SP_Matrix_Data * D =  (SP_Matrix_Data *)calloc( 1, sizeof(SP_Matrix_Data) );
  int j;
  
  D->Presence = (double **)calloc(No_of_SITES, sizeof(double *));
  for (j=0; j<No_of_SITES; j++)
      D->Presence[j] = (double *)calloc(No_of_TIMES, sizeof(double ));
  
  D->Time_Vector = (double *)calloc(Total_No_of_TRANSECTS, sizeof(double));
  D->Transects   = (int *)calloc(No_of_TIMES, sizeof(int));

  D->Sp_Time = (double **)calloc(No_of_SITES, sizeof(double *));
  for (j=0; j<No_of_SITES; j++)
      D->Sp_Time[j] = (double *)calloc( No_of_TIMES, sizeof(double) );
  
  D->No_Sp_Time = (int *)calloc( No_of_SITES, sizeof(int) );

  D->Sp_Transects = (int **)calloc(No_of_SITES, sizeof(int *)); 
  for (j=0; j<No_of_SITES; j++)
    D->Sp_Transects[j]  = (int *)calloc( No_of_TIMES, sizeof(int) );

  D->Sp_Total_No_Transects = (int *)calloc(No_of_SITES, sizeof(int *)); 
  
  D->Name = (char *)calloc(50, sizeof(char) );
  
  return(D);
}				       

void SP_Matrix_Data_Free( SP_Matrix_Data * D)
{
  int j; 
  int No_of_SITES; 
 
  No_of_SITES = D->No_of_SITES;
  
  for (j=0; j<No_of_SITES; j++) free( D->Presence[j] );
  free( D->Presence );
  
  for (j=0; j<No_of_SITES; j++) free( D->Sp_Time[j] );
  free( D->Sp_Time );

  for (j=0; j<No_of_SITES; j++) free( D->Sp_Transects[j] );
  free( D->Sp_Transects );
  
  free (D->Sp_Total_No_Transects);
  free (D->Time_Vector);
  free (D->Transects);
  free (D->No_Sp_Time);
  free (D->Name);   
  free (D);
}
  
void SP_Matrix_Data_Setup(int No_of_SITES, int No_of_COLUMNS,
			  int Total_No_of_TRANSECTS, 
			  SP_Matrix_Data * D, double ** Presence, 
			  double * Time_Vector, double ** Sp_Time, int * No_Sp_Time, 
			  int * Transects, 
			  char * Sp_Name )
{
  int i,j;

  D->No_of_SPECIES = No_of_SITES; // Species: Number of rows of the presence atrix.
  D->No_of_SITES = No_of_SITES;   // Sites: Number of rows of the presence matrix
                                  // No_of_SITES and No_of_SPECIES are synonimous 
                                  // interchangeable members of the data structure
                                  // Their values should match.
  D->No_of_TIMES = No_of_COLUMNS;           
  // No of sampling times: Number of columns (if there are no replicates per 
  // sampling time).
  D->Total_No_of_TRANSECTS = Total_No_of_TRANSECTS; 
  // This variable matches No_of_TIMES if there are no replicates per sampling time 

  for (i=0; i<No_of_SITES; i++){
    D->No_Sp_Time[i] = No_Sp_Time[i];
    for (j=0; j<No_of_COLUMNS; j++) {
      D->Presence[i][j] = Presence[i][j];
      D->Sp_Time[i][j]  = Sp_Time[i][j]; 
    }
  }

  for (j=0; j<No_of_COLUMNS; j++) D->Transects[j] = Transects[j]; 
  
  for (j=0; j<Total_No_of_TRANSECTS; j++) D->Time_Vector[j] = Time_Vector[j]; 

  memcpy( D->Name, Sp_Name, strlen(Sp_Name)+1 );
}

void SP_Matrix_Data_Uneven_Setup( SP_Matrix_Data * D, double ** Presence, 
				 int No_of_SITES, int Total_No_of_TRANSECTS, 
				 double * Time_Vector, int * Transects,
				 int No_of_COLUMNS,
				 double ** Sp_Time, int * No_Sp_Time,
				 int ** Sp_Transect, int * Sp_No_of_Transects )
{
  int i,j;

  D->No_of_SPECIES = No_of_SITES; // Species: Number of rows of the presence atrix.
  D->No_of_SITES = No_of_SITES;   // Sites: Number of rows of the presence matrix
                                  // No_of_SITES and No_of_SPECIES are synonimous 
                                  // interchangeable members of the data structure
                                  // Their values should match.
  D->No_of_TIMES = No_of_COLUMNS;           
  // No of sampling times: Number of columns (if there are no replicates per 
  // sampling time).
  D->Total_No_of_TRANSECTS = Total_No_of_TRANSECTS; 
  // This variable matches No_of_TIMES if there are no replicates per sampling time 

  for (i=0; i<No_of_SITES; i++){
    D->No_Sp_Time[i] = No_Sp_Time[i];
    D->Sp_Total_No_Transects[i] = Sp_No_of_Transects[i]; 
    for (j=0; j<Sp_No_of_Transects[i]; j++) {
      D->Presence[i][j] = Presence[i][j];
    }
    for (j=0; j<No_Sp_Time[i]; j++) {
      D->Sp_Time[i][j]  = Sp_Time[i][j];
      D->Sp_Transects[i][j] = Sp_Transect[i][j]; 
    }
  }

  for (j=0; j<No_of_COLUMNS; j++) {
    D->Transects[j] = Transects[j];
    D->Time_Vector[j] = Time_Vector[j]; 
  }

}
  
void SP_Matrix_Data_Write( SP_Matrix_Data * D )
{
  /* This function does not works at all unless the 'D' 
     data structure has been fully set up */
  
  int j,k,n,m;
  
  int No_of_SITES           = D->No_of_SITES;
  // int Total_No_of_TRANSECTS = D->Total_No_of_TRANSECTS;
  // int No_of_TIMES           = D->No_of_TIMES;

  for(j=0; j<No_of_SITES; j++) {
    
    printf(" Sampling Times (%d-th row) = %d\t Time(No of Transects) = {",
	   j, D->No_Sp_Time[j] );
    for(k=0; k<D->No_Sp_Time[j]; k++)
      printf(" %g(%d) ", D->Sp_Time[j][k], D->Sp_Transects[j][k]);
    printf("}\n");
    
    m = 0;
    printf(" { "); 
    for(k=0; k<D->No_Sp_Time[j]; k++)
      for(n=0; n<D->Sp_Transects[j][k]; n++) printf("%g ", D->Presence[j][m++]);
    printf("}\n"); 
  }
  printf("\n\n");
}
