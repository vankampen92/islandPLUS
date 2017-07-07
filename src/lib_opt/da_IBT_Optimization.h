/* Model Seletion Functions */

void MODEL_SELECTION_UPGMA_R_SHLIB( double Colonization_Rate, double * C_Range, 
				    double Extinction_Rate, double * E_Range, 
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
				    double ** Results );

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
					     double ** Results ); 

/* Functions that use colonization and extinction probabilities 
   rather than rates are only possible when sampling 
   times are distributed at regular intervals at least for
   a good fraction of transitions */

int Checking_for_Parameter_Correctness( double Colonization_Rate, 
					double Extinction_Rate,
					double Detectability_Value, 
					double Phi_0 );

double NLLikelihood_Calculation ( int n, Time_Control * Time, 
				  double ** Presence_Data, int No_of_SPECIES,
				  double Colonization, double Extinction );

double NLLikelihood_Calculation_Transition_Probabilities ( int n, Time_Control * Time, 
							   double ** Presence_Data, 
							   int No_of_SPECIES,
							   double Colonization_Probability, 
							   double Extinction_Probability );

/* Functions using colonization and extinction temporal rates 
   rather than probabilities. These functions do not have any 
   restriction (unlike above) so sampling times can be distributed 
   in any possible way */

/* Functions using the 2-parameter model (extinction and colonization) */
double GSL_NLLikelihood_Function ( const gsl_vector * x, void * Par );

void mle_Col_Ext_Uneven_Matrix_R_SHLIB (double ** Presence, int S, 
					double * Time_Vector, int No_of_TIMES, 
					double ** Sp_Time_Vector, int * T, 
					double * Colonization, double * C_Range,
					double * Extinction,  double * E_Range, 
					int * No_of_PARAMETERS,
					int * No_of_PARAMETERS_MAX,
					int * Index, 
					int * Discretization, 
					double * Tolerance, 
					int * No_of_ITERATIONS,
					int * Verbose, 
					int * Minimization,
					double * Value );

double GSL_NLL_Function_Uneven( const gsl_vector * x, void * Par );

/* Functions using the 4-parameter model (extinction, colonization, detectability, 
	                                  phi_0) */
double NLL_MacKenzie_000_000(double Psi, double a, double C, double E ); 

double GSL_MacKenzie_NLLikelihood_Function ( const gsl_vector * x, void * Par );

double GSL_MacKenzie_NLL_Uneven_Function ( const gsl_vector * x, void * Par ); 

double MacKenzie_NLLikelihood_Calculation ( double ** Presence_Data, int No_of_SPECIES, int N, 
					    double * Time_Vector, int * Transects, int n, 
					    double Colonization_Rate, 
					    double Extinction_Rate,
					    double Detectability_Value, 
					    double Phi_0 );

/* Functions that are used for both models */
double GSL_Minimization_Simplex (Parameter_Fitting * F, 
				 gsl_vector * Initial_Guess, 
				 gsl_vector * Solution, 
				 double ( * Function)( const gsl_vector * , void * ) );

void mle_MacKenzie_Uneven_Matrix_R_SHLIB( double ** Presence, int S, int N, 
					  double * Time_Vector, int * Transects, int T, 
					  double ** Sp_Time_Vector, int * Sp_T, 
					  int ** Sp_Transects, int * Sp_N, 
					  double * Colonization, double * C_Range, 
					  double * Extinction, double * E_Range, 
					  double * Detectability, double * D_Range, 
					  double * Phi_Time_0, double * P_Range,
					  int * No_of_PARAMETERS,
					  int * No_of_PARAMETERS_MAX, 
					  int * Index, int * Discretization,
					  double * Tolerance,
					  int * No_of_ITERATIONS,
					  int * Verbose, 
					  int * Minimization,
					  double * Value ); 


