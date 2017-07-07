#include "HEADERS.h"
#include "da_IBT_Functions.h"

void AIC_Summary_RESULTS( int * No_of_MODEL_PARAMETERS, int No_of_SPECIES, 
			  double * MODEL_NLL, double * MODEL_AIC, double * MODEL_AIC_c, 
			  double ** Results )
{ 
  int i;
  double * AIC_d         = (double *)calloc( No_of_SPECIES, sizeof(double) );
  double * AIC_w         = (double *)calloc( No_of_SPECIES, sizeof(double) );
  double D, E;
  
  double min_AIC_c = MODEL_AIC_c[0];
  int min_i        = 0; 
  for(i=0; i<No_of_SPECIES; i++) {
    if( MODEL_AIC_c[i] < min_AIC_c) min_i = i;
    min_AIC_c = MIN( MODEL_AIC_c[i], min_AIC_c );
  }
  
  E = 0.0; 
  for(i=0; i<No_of_SPECIES; i++) { 
    D = AIC_d[i] = MODEL_AIC_c[i] - min_AIC_c;
    E += exp( - 0.5 * D ); 
  }
  for(i=0; i<No_of_SPECIES; i++){  
    D = AIC_d[i]; 
    AIC_w[i] = exp( - 0.5 * D ) / E; 
  }
  
  for(i=0; i<No_of_SPECIES; i++) {  
    Results[i][0] = (double)No_of_MODEL_PARAMETERS[i]; 
    Results[i][1] = MODEL_NLL[i];	
    Results[i][2] = MODEL_AIC[i];
    Results[i][3] = MODEL_AIC_c[i];
    Results[i][4] = AIC_d[i];
    Results[i][5] = AIC_w[i];
  }

  free(AIC_d);
  free(AIC_w);
}

void Model_Selection_AIC_Latex_Table( char ** Name, 
				      int *** PARTITION, int * G, int ** K, 
				      int No_of_SPECIES, 
				      double * NLL, double ** COL, double ** EXT, 
				      double * AIC, double * AIC_c )
{
  int i,j,k;
  double * AIC_d         = (double *)calloc( No_of_SPECIES, sizeof(double) );
  double * AIC_w         = (double *)calloc( No_of_SPECIES, sizeof(double) );
  double D, E;
  char *  p;
  
  double min_AIC_c = AIC_c[0];
  int min_i        = 0; 
  for(i=0; i<No_of_SPECIES; i++) {
    if( AIC_c[i] < min_AIC_c) min_i = i;
    min_AIC_c = MIN( AIC_c[i], min_AIC_c );
  }
 
  E = 0.0; 
  for(i=0; i<No_of_SPECIES; i++) { 
    D = AIC_d[i] = AIC_c[i] - min_AIC_c;
    E += exp( - 0.5 * D ); 
  }
  for(i=0; i<No_of_SPECIES; i++){  
    D = AIC_d[i]; 
    AIC_w[i] = exp( - 0.5 * D ) / E; 
  }

  for(i=0; i<No_of_SPECIES; i++) {
    Rprintf(" Partition %d-th: Number of estimated parameters: %d\n", i, G[i]*2);
    for(j=0; j<G[i]; j++) { 
	Rprintf("{ ");
	for(k=0; k<K[i][j]; k++) Rprintf("%s ", Name[PARTITION[i][j][k]]);
	Rprintf("} ");
    }
    Rprintf("\n");
    Rprintf(" NLL = %g\t AIC = %g\tAIC (corrected) = %g\t", NLL[i], AIC[i], AIC_c[i]); 	  
    Rprintf(" AIC_d = %g\t AIC_w = %g\n", AIC_d[i], AIC_w[i] );
  }

  double ** VALUE = (double **)calloc( No_of_SPECIES, sizeof(double *) );
  for(i=0; i<No_of_SPECIES; i++) { 
    VALUE[i] = (double *)calloc( 5, sizeof(double) );
    
    VALUE[i][0] = NLL[i];
	
    VALUE[i][1] = AIC[i];

    VALUE[i][2] = AIC_c[i];

    VALUE[i][3] = AIC_d[i];

    VALUE[i][4] = AIC_w[i];
  }
  char * Num = (char *)calloc( 10, sizeof(char) );
  char ** Row_Name = (char **)calloc( No_of_SPECIES, sizeof(char *) );
  int No_of_PARAMETERS;
  for(i=0; i<No_of_SPECIES; i++) { 
    Row_Name[i] = (char *)calloc( 20, sizeof(char) );
    Row_Name[i][0] = '\0';
    No_of_PARAMETERS = 2*G[i]; 
    sprintf( Num, "%d", No_of_PARAMETERS );
    p = strcat( Row_Name[i], Num );
    p = strcat( Row_Name[i], "-parameter model" ); 
  }
  char ** Column_Name = (char **)calloc( 6, sizeof(char *) );
  for(i=0; i<6; i++) {
    Column_Name[i] = (char *)calloc( 20, sizeof(char) );
    Column_Name[i][0] = '\0';
    switch (i) {
    case 0: p = strcat( Column_Name[i], "Model" );
      break;
    case 1: p = strcat( Column_Name[i], "NLL" );
      break;
    case 2: p = strcat( Column_Name[i], "AIC" );
      break;
    case 3: p = strcat( Column_Name[i], "AIC corrected" );
      break;
    case 4: p = strcat( Column_Name[i], "AIC difference" );
      break;
    case 5: p = strcat( Column_Name[i], "AIC weights" );
      break;
    default: Rprintf(" Index j = %d out of range (0,...,5) (Model Selection Latex Table Function)\n", j);
      error(0,0,"Program aborted");
    }
  }
  Latex_Table_Driver( "Model_Selection_Results.tex", 
		      No_of_SPECIES, 6, Row_Name, Column_Name, VALUE );
  
  for(i=0; i<6; i++) free(Column_Name[i]);
  free(Column_Name);
  for(i=0; i<No_of_SPECIES; i++) { free(VALUE[i]); free(Row_Name[i]); }
  free(VALUE); free(Row_Name);
  free(Num);

  /* Colonization - Extinction parameters of the best model */
  VALUE = (double **)calloc( G[min_i], sizeof(double *) );
  for(j=0; j<G[min_i]; j++) { 
    VALUE[j] = (double *)calloc( 2, sizeof(double) );
    
    VALUE[j][0] = EXT[min_i][j];
    VALUE[j][1] = COL[min_i][j];
  }
  Column_Name = (char **)calloc( 3, sizeof(char *) );
  for(i=0; i<3; i++) {
    Column_Name[i] = (char *)calloc( 20, sizeof(char) );
    Column_Name[i][0] = '\0';
    switch (i) {
    case 0: p = strcat( Column_Name[i], "Species Group" );
      break;
    case 1: p = strcat( Column_Name[i], "Extinction Rate" );
      break;
    case 2: p = strcat( Column_Name[i], "Colonization Rate" );
      break;
    default: Rprintf(" Index i = %d out of range (0,1,2) (Model Selection Latex Table Function)\n", j);
      error(0,0,"Program aborted");
    }
  }
  Row_Name = (char **)calloc( G[min_i], sizeof(char *) );
  for(j=0; j<G[min_i]; j++) { 
    Row_Name[j] = (char *)calloc( 50, sizeof(char) );
    Row_Name[j][0] = '\0';
    p = strcat( Row_Name[j], "{ ");
    for(k=0; k<K[min_i][j]; k++) { 
      p = strcat( Row_Name[j], Name[PARTITION[min_i][j][k]]);
      p = strcat( Row_Name[j], " ");
    }
    p = strcat( Row_Name[j], " }" );
  }
  Latex_Table_Driver ( "Best_Model_Colonization_Extinction_Results.tex", 
		       G[min_i], 3, Row_Name, Column_Name, VALUE );
  
  for(i=0; i<3; i++) free(Column_Name[i]);
  free(Column_Name);
  for(i=0; i<G[min_i]; i++) { free(VALUE[i]); free(Row_Name[i]); }
  free(VALUE); free(Row_Name);

  free(AIC_d);
  free(AIC_w);
}

void Latex_Table_Driver (char * Name_of_File, 
			 int No_of_ROWS, int No_of_COLUMNS, 
			 char ** Row_Name, char ** Column_Name, 
			 double ** VALUE )
{
  int i,j;
  FILE * fp;
 
  fp = fopen(Name_of_File, "w");
  fprintf(fp, "\\input{TableOpening}\n");

  fprintf(fp, "\\begin{table}\n");
  fprintf(fp, "   \\centering\n");
  fprintf(fp, "   \\begin{tabular}{l");
  for(i=1; i<No_of_COLUMNS; i++) fprintf(fp, "c");
  fprintf(fp, "}\n");
  fprintf(fp, "%s", Column_Name[0]);
  for(i=1; i<No_of_COLUMNS; i++) 
    fprintf(fp, "& %s", Column_Name[i]);
  fprintf(fp, "\\"); fprintf(fp, "\\"); fprintf(fp, "\n");
  fprintf(fp, "\\hline\n"); 
  for(i=0; i<No_of_ROWS; i++) {
    fprintf(fp, "%s", Row_Name[i]);
    for(j=1; j<No_of_COLUMNS; j++) 
      fprintf(fp, "& %g", VALUE[i][j-1]);
    fprintf(fp, "\\"); fprintf(fp, "\\"); fprintf(fp, "\n");
  }
  fprintf(fp, "   \\end{tabular}\n");
  fprintf(fp, "   \\caption{Caption goes here}\n");
  fprintf(fp, "   \\label{tab:myfirsttable}\n");
  fprintf(fp, "\\end{table}\n");

  fprintf(fp, "\\end{document}\n");
  fclose(fp);
}
