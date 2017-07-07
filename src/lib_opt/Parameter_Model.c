#include <MODEL.h>

// R_SHLIB is a compilation FLAG for the generation of libraries compatible 
// with the generation of a shared object to be a callable R function. 
// When no shared object is needed, the R_SHLIB flag is not required.

#ifndef R_SHLIB

#include <include.Parameter_Model.extern.h>
void  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N ( Parameter_Model * P )
{
  P->No_of_SPECIES     = No_of_SPECIES;
  
  P->No_of_COLUMNS     = No_of_COLUMNS;

  P->No_of_TRANSECTS   = No_of_TRANSECTS;
  
  P->Colonization_Rate = Colonization_Rate;     /* Key 0 */

  P->Extinction_Rate   = Extinction_Rate;       /* Key 1 */

  P->Detectability_Value = Detectability_Value; /* Key 2 */

  P->Phi_0             = Phi_0;                 /* Key 3 */

  P->RATES             = RATES;
}
#else

void Parameter_Model_Initialization_R_SHLIB ( Parameter_Model * P,
					      int S, int T,
					      double * Colonization,
					      double * Extinction, 
					      double * Detectability,
					      double * Phi_Time_0 )
{
  P->No_of_SPECIES     = S;
  
  P->No_of_COLUMNS     = T;

  P->Colonization_Rate = (* Colonization );     /* Key 0 */

  P->Extinction_Rate     = (* Extinction );     /* Key 1 */

  P->Detectability_Value = (* Detectability );  /* Key 2 */

  P->Phi_0   = (* Phi_Time_0 );                 /* Key 3 */
 
  P->RATES   = 1;  /* Extinction and Colonization should be true rates */
}
#endif

void P_A_R_A_M_E_T_E_R___M_O_D_E_L___F_R_E_E( Parameter_Model * P)
{
  free ( P );
}

void Parameter_Model_Initialization_From_Values ( Parameter_Model * P, 
						  double pColonization_Rate, 
						  double pExtinction_Rate, 
						  double pDetectability_Value, 
						  double pPhi_0 )
{
  P->Colonization_Rate = pColonization_Rate;      /* Key 0 */

  P->Extinction_Rate   = pExtinction_Rate;        /* Key 1 */

  P->Detectability_Value = pDetectability_Value;  /* Key 2 */

  P->Phi_0   = pPhi_0;                            /* Key 3 */
}

void Vector_Entries_into_Parameter_Model ( const gsl_vector * X, Parameter_Model * P, 
					   int * Parameter_Index, int No_of_PARAMETERS )
{
  int i; 
  int key;
  double value; 

  for( i=0; i<No_of_PARAMETERS; i++) { 
    key = Parameter_Index[i];
    value = gsl_vector_get(X, i);
    Vector_Entry_into_Parameter_Model ( value, key, P );
  }
}

void Vector_Entry_into_Parameter_Model ( double value, 
					 int key, Parameter_Model * P )
{
  switch (key) {
  case 0:
    P->Colonization_Rate = value;
    break;
  case 1:
    P->Extinction_Rate = value;
    break;
  case 2:
    P->Detectability_Value = value;
    break;
  case 3:
    P->Phi_0 = value;
    break;
  default:
    Rprintf(" Error from:\n"); 
    Rprintf(" Vector_Entry_into_Parameter_Model(...) in Parameter_Model.c\n");
    Rprintf(" INVALID PARAMETER KEY [key = %d]\n", key);
    int N = MODEL_PARAMETERS_MAXIMUM;
    Rprintf(" The maximum number of parameters is %d\n", N);
    Rprintf(" The permited keys go from 0, to %d\n", N-1);
    error(0,0,"Program aborted");
  }
}

void Parameter_Model_into_Vector_Entries ( Parameter_Model * P, gsl_vector * X,  
					   int * Parameter_Index, int No_of_PARAMETERS )
{
  int i; 
  int key;
  double value; 

  for( i=0; i<No_of_PARAMETERS; i++) { 
    key = Parameter_Index[i];    
    value = Parameter_Model_into_Vector_Entry( key, P );
    gsl_vector_set(X, i, value);
  }
}

double Parameter_Model_into_Vector_Entry ( int key, Parameter_Model * P )
{
  double value;

  switch (key) {
  case 0:
    value = P->Colonization_Rate;
    break;
  case 1:
    value = P->Extinction_Rate;
    break;
  case 2:
    value = P->Detectability_Value;
    break;
  case 3:
    value = P->Phi_0;
    break;
  default:
    Rprintf(" Error from:\n"); 
    Rprintf(" Parameter_Model_into_Vector_Entry(...) in Parameter_Model.c\n");
    Rprintf(" INVALID PARAMETER KEY [key = %d]\n", key);
    int N = MODEL_PARAMETERS_MAXIMUM;
    Rprintf(" The maximum number of parameters is %d\n", N);
    Rprintf(" The permited keys go from 0, to %d\n", N-1);
    error(0,0,"Program aborted");
  }

  return(value);
}
