void Latex_Table_Driver (char * Name_of_File, 
			 int No_of_ROWS, int No_of_COLUMNS, 
			 char ** Row_Name, char ** Column_Name, 
			 double ** VALUE );

void Model_Selection_AIC_Latex_Table( char ** Name, 
				      int *** PARTITION, int * G, int ** K, 
				      int No_of_SPECIES, 
				      double * NLL, double ** COL, double ** EXT, 
				      double * AIC, double * AIC_c ); 

void AIC_Summary_RESULTS( int * No_of_MODEL_PARAMETERS, int No_of_SPECIES,  
			  double * MODEL_NLL, double * MODEL_AIC, double * MODEL_AIC_c, 
			  double ** Results);

