#include <MODEL.h>

// R_SHLIB is a compilation FLAG for the generation of libraries compatible 
// with the generation of a shared object. When no shared object is needed, 
// the R_SHLIB flag is not required.

#if defined R_SHLIB
void Parameter_Space_Alloc_R_SHLIB( Parameter_Space * S, int No_of_PAR, 
				    int No_0, int No_1, int No_2, int No_3 )
{
  S->P_MAX    = gsl_vector_alloc(No_of_PAR);
  S->P_min    = gsl_vector_alloc(No_of_PAR);
  S->Accuracy = gsl_vector_alloc(No_of_PAR);

  S->Parameter_Index =(int *)calloc(No_of_PAR, sizeof(int) );
  
  S->Parameter_MAX = (Parameter_Model *)calloc(1, sizeof(Parameter_Model) ) ;
  S->Parameter_min = (Parameter_Model *)calloc(1, sizeof(Parameter_Model) ) ;
  S->Parameter_Acc = (Parameter_Model *)calloc(1, sizeof(Parameter_Model) ) ;

  S->N = (int *)calloc(No_of_PAR, sizeof(int) );
  S->N_Par_Value = (double **)calloc(No_of_PAR, sizeof(double *) );
  
  if( No_of_PAR > 0 ) S->N_Par_Value[0] = (double *)calloc( No_0, sizeof(double) ); 
  if( No_of_PAR > 1 ) S->N_Par_Value[1] = (double *)calloc( No_1, sizeof(double) );
  if( No_of_PAR > 2 ) S->N_Par_Value[2] = (double *)calloc( No_2, sizeof(double) );
  if( No_of_PAR > 3 ) S->N_Par_Value[3] = (double *)calloc( No_3, sizeof(double) ); 
}

void Parameter_Space_Boundaries_R_SHLIB( Parameter_Space * S, 
					 double * Range_COL, 
					 double * Range_EXT,
 					 double * Range_DET,
	                                 double * Range_PHI )
{
  S->Parameter_MAX->Colonization_Rate   = Range_COL[1];   /* Key 0 */
  S->Parameter_MAX->Extinction_Rate     = Range_EXT[1];   /* Key 1 */
  S->Parameter_MAX->Detectability_Value = Range_DET[1];   /* Key 2 */
  S->Parameter_MAX->Phi_0               = Range_PHI[1];   /* Key 3 */
  
  S->Parameter_min->Colonization_Rate   = Range_COL[0];   /* Key 0 */
  S->Parameter_min->Extinction_Rate     = Range_EXT[0];   /* Key 1 */
  S->Parameter_min->Detectability_Value = Range_DET[0];   /* Key 2 */
  S->Parameter_min->Phi_0               = Range_PHI[0];   /* Key 3 */
}

void Parameter_Space_Accuracies_R_SHLIB( Parameter_Space * S, 
					 double Acc_C, double Acc_E, double Acc_D, double Acc_P  )
{
  Parameter_Model_Initialization_From_Values(S->Parameter_Acc, 
					     Acc_C, Acc_E, Acc_D, Acc_P );
}

void Parameter_Space_Initialization_R_SHLIB(Parameter_Space * S, 
					    double TOL, 
					    int No_of_ITERATIONS,
					    int No_of_PARAMETERS, 
	 				    int I0, int I1, int I2, int I3,
					    int No_0, int No_1, int No_2, int No_3)
{
  int i, j, key;
  double lo_P, hi_P, ac_P;

  S->TOLERANCE        = TOL;
  S->MAX_No_of_ITERATIONS = No_of_ITERATIONS;
  S->No_of_PARAMETERS =  No_of_PARAMETERS;

  if( No_of_PARAMETERS > 0 ) S->Parameter_Index[0] = I0; S->N[0] = No_0; 
  if( No_of_PARAMETERS > 1 ) S->Parameter_Index[1] = I1; S->N[1] = No_1;
  if( No_of_PARAMETERS > 2 ) S->Parameter_Index[2] = I2; S->N[2] = No_2; 
  if( No_of_PARAMETERS > 3 ) S->Parameter_Index[3] = I3; S->N[3] = No_3; 
  
  for ( i=0; i<S->No_of_PARAMETERS; i++ ) {  
      key = S->Parameter_Index[i];
      lo_P = Parameter_Model_into_Vector_Entry( key, S->Parameter_min );
      hi_P = Parameter_Model_into_Vector_Entry( key, S->Parameter_MAX );	
      ac_P = Parameter_Model_into_Vector_Entry( key, S->Parameter_Acc );
      gsl_vector_set(S->P_min, i, lo_P); 
      gsl_vector_set(S->P_MAX, i, hi_P);
      gsl_vector_set(S->Accuracy, i, ac_P);
      
      for(j=0; j<S->N[i]; j++){
	S->N_Par_Value[i][j] = lo_P + (double)j * (hi_P - lo_P)/(double)(S->N[i]-1);
      }
  }                                                                           
}

void Parameter_Space_Free(Parameter_Space * S)
{
  gsl_vector_free(S->P_MAX);
  gsl_vector_free(S->P_min);
  gsl_vector_free(S->Accuracy);
  free (S->Parameter_Index);
  free (S->Parameter_MAX);
  free (S->Parameter_min);
  free (S->Parameter_Acc);

  if( S->No_of_PARAMETERS > 0 ) free( S->N_Par_Value[0] );
  if( S->No_of_PARAMETERS > 1 ) free( S->N_Par_Value[1] );
  if( S->No_of_PARAMETERS > 2 ) free( S->N_Par_Value[2] );
  if( S->No_of_PARAMETERS > 3 ) free( S->N_Par_Value[3] );
  
  free( S->N_Par_Value );
  free( S->N );

  free( S );
}
#endif
