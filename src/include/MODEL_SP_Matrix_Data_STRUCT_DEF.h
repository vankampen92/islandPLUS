typedef struct SP_Matrix_Datainfo
{
  /* * * Model Parameters  * * */
  double ** Presence;
  int No_of_SPECIES;         // Species: Number of rows of the presence atrix.
  int No_of_SITES;           // Sites: Number of rows of the presence matrix
                             // No_of_SITES and No_of_SPECIES are synonymous 
                             // interchangeable members of the data structure
                             // Their values should match.
  int No_of_TIMES;           // No of sampling times: Number of columns (if there are
                             // no replicates per sampling time).
  int Total_No_of_TRANSECTS; // This variable matches No_of_TIMES if there are
                             // no replicates per sampling time 
  double * Time_Vector;
  int    * Transects;

  double ** Sp_Time;         // This is to deal with different times sampling schemes 
  int    ** Sp_Transects;    // per row both times and number of transects (in case  
  int    *  No_Sp_Time;      // there are more than 1 transect per sampling time)
  int    *  Sp_Total_No_Transects; 
  int   *** Structure ;      // an 1/0 recording the missing value structure of 
                             // the initial data matrix. 
 
  char * Name;              // Name of the species/guild/Phylum stored in this 
                            // structure      
}SP_Matrix_Data;


/* The following functions handle this data structure and are compiled in the library:
// libda_IBT_Optimization.a (see ./Library_Optimization/SP_Matrix_Data.c)
*/

SP_Matrix_Data * SP_Matrix_Data_Alloc( int , int , int ); 

void SP_Matrix_Data_Setup(int , int , int , 
			 SP_Matrix_Data * , double ** , 
			 double * , double ** , int * , int * , 
			 char * ); 

void SP_Matrix_Data_Uneven_Setup( SP_Matrix_Data * , double ** , 
				 int , int , 
				 double * , int * ,
				 int ,
				 double ** , int * ,
				 int ** , int * ) ;

void SP_Matrix_Data_Free( SP_Matrix_Data * );

void SP_Matrix_Data_Write( SP_Matrix_Data * );
