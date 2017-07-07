typedef struct Parameter_Spaceinfo
{
  /* Model Parameters defining the boundary of 
     the Parameter Space */
  Parameter_Model * Parameter_min;
  Parameter_Model * Parameter_MAX;
  
  /* Model Parameters defining parameter dependent accuracies */
  Parameter_Model * Parameter_Acc;

  /* * * Model Parameters for local community dynamics  * * */
  gsl_vector * P_MAX;         /* Max and Min for each parameter */
  gsl_vector * P_min;     
  gsl_vector * Accuracy;      /* Different accuracy for each searcheable parameter */
  
  
  int * Parameter_Index;

  /* This will define GRID on which different parameter combinations will be 
     evalutated */
  int * N; 
  double ** N_Par_Value;

  #include <include.Parameter_Space.global.h>
 
}Parameter_Space;

// R_SHLIB is a compilation FLAG for the generation of libraries compatible 
// with the generation of a shared object. When no shared object is needed, 
// the R_SHLIB flag is not required.
// Functions at work when R_SHLIB is NOT defined 

// These functions are compiled in the library:
// libda_IBT_Optimization.a
// (see ./lib_opt/Parameter_Space.c)

void Parameter_Space_Alloc( Parameter_Space * S );
void Parameter_Space_Initialization(Parameter_Space * S);
void Parameter_Space_Boundaries( Parameter_Space * S );
void Parameter_Space_Accuracies( Parameter_Space * S );
// Functions at work when R_SHLIB is defined 
void Parameter_Space_Alloc_R_SHLIB(Parameter_Space * , int , int , int , int , int  );
void Parameter_Space_Boundaries_R_SHLIB(Parameter_Space * , double * , double * , double * , double * );
void Parameter_Space_Accuracies_R_SHLIB( Parameter_Space *, double , double , double , double );
void Parameter_Space_Initialization_R_SHLIB(Parameter_Space * , 
					    double , 
					    int ,
					    int , 
					    int , int, int , int,
					    int , int, int , int );
// Functions at work either way
void Parameter_Space_Free(Parameter_Space * );
void Parameter_Index_Checking_Ordering (int * , int * , int , 
					int * , int * , int * , int * );

