void IO_Filtering_Out_Missing_Values ( int No_of_SPECIES,  
				       double *** Presence, int * No_of_SITES, 
				       double ** Time_Vector, int * No_of_TIMES, 
				       double *** Sp_Time, int ** No_of_Sp_Times , 
				       double MISSING_VALUE_FLAG );

void IO_Filtering_Out_Matrix( double ** Presence, int * No_of_SPECIES,  
			      double * Time, int No_of_TIMES,    
			      double ** Sp_Time, int * No_of_Sp_Times, 
			      double MISSING_VALUE_FLAG );

