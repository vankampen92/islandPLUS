typedef struct Parameter_Modelinfo
{
  /* * * Model Parameters  * * */
#include <include.Parameter_Model.global.h>  

  Time_Control * Time;

}Parameter_Model;

// R_SHLIB is a compilation FLAG for the generation of libraries compatible 
// with the generation of a shared object. When no shared object is needed, 
// the R_SHLIB flag is not required.

// These functions are compiled in the library:
// libda_IBT_Optimization.a
// (see ./lib_opt/Parameter_Model.c)

void  P_A_R_A_M_E_T_E_R___I_N_I_T_I_A_L_I_Z_A_T_I_O_N (Parameter_Model * P);
void Parameter_Model_Initialization_R_SHLIB ( Parameter_Model *, 
					      int , int , 
					      double * , double *, double *, double * );

void P_A_R_A_M_E_T_E_R___M_O_D_E_L___F_R_E_E( Parameter_Model * P);

void Parameter_Model_Initialization_From_Values ( Parameter_Model * , 
						  double , 
						  double , 
						  double ,
						  double );

void Vector_Entries_into_Parameter_Model ( const gsl_vector * X, Parameter_Model * P, 
					   int * Parameter_Index, int No_of_PARAMETERS );

void Vector_Entry_into_Parameter_Model ( double value, 
					 int key, Parameter_Model * P );

void Parameter_Model_into_Vector_Entries ( Parameter_Model * P, gsl_vector * X,  
					   int * Parameter_Index, int No_of_PARAMETERS );

double Parameter_Model_into_Vector_Entry ( int key, Parameter_Model * P );
