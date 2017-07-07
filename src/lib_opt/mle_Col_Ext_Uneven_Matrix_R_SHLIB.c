#include <MODEL.h>

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
					double * Value )
{
  /* This function is just a wrapper for a GSL minimization algorithm:
 
        GSL_Minimization_Simplex( Parameter_Fitting * F, 
				  gsl_vector * Initial_Guess, 
				  gsl_vector * Solution )
				  
     The input matrix, however, is uneven, this is, each 
     presence/absence row has a different length and a 
     different, characteristic sampling vector. The Presence[][] 
     has, then, the folowing structure:

      Sp
      ---
      1:  0  1  0  0  0  1  (associated to vector Sp_Time_Vector[0])   
      2:  0  1 
      ...
      S:  1  1  0  1        (associated to vector Sp_Time_vector[S-1])

     Input arguments:
     .  T[], number of temporal observations per species (no of rows)
     .  S, number of species (no of rows in Presence matrix)
     .  Presence, presence 1/0 data matrix. 
     .  Sp_Time_Vector[], double array containing sampling times per
     species (for each row).
     
     Output arguments:
     .  Colonization, maximum likelihood estimate of the colonization rate
     .  Extinction, maximum likelihood estimate of the extinction rate
     .  Value, the value of the NLL at the point of the Colonization and Extinctin 
     mle estimates
  */
  SP_Matrix_Data  * Data = (SP_Matrix_Data *)calloc( 1, sizeof(SP_Matrix_Data) );
  Data->Presence = Presence;
  Data->No_of_SPECIES = S;
  Data->No_of_TIMES   = No_of_TIMES;
  Data->Time_Vector   = Time_Vector;
  Data->Sp_Time    = Sp_Time_Vector;
  Data->No_Sp_Time = T;

  Parameter_Model * P = (Parameter_Model *)calloc(1, sizeof(Parameter_Model) );
  P->No_of_SPECIES       = S;
  P->No_of_COLUMNS       = No_of_TIMES;
  P->No_of_TRANSECTS     = No_of_TIMES;
  P->Colonization_Rate   = (* Colonization);      /* Key 0 */
  P->Extinction_Rate     = (* Extinction);        /* Key 1 */
  P->Detectability_Value = 1.0;                   /* Key 2 */
  P->Phi_0               = 1.0;                   /* Key 3 */
  P->RATES               = 1;
  
  Parameter_Space * Space = (Parameter_Space *)calloc(1, sizeof(Parameter_Space) );
  Parameter_Space_Alloc_R_SHLIB( Space, (* No_of_PARAMETERS), 
				 Discretization[0], Discretization[1], 
				 Discretization[2], Discretization[3] );

  if( (* No_of_PARAMETERS_MAX) != (* No_of_PARAMETERS) ) error(0,0,"Number of Parameters do not match"); 
  int No_C, No_E, No_D, No_P; 
  Parameter_Index_Checking_Ordering(Index, Discretization, (* No_of_PARAMETERS_MAX), 
				    &No_C, &No_E, &No_D, &No_P );
  double Acc_C, Acc_E, Acc_D, Acc_P;
  Acc_C = (C_Range[1] - C_Range[0])/((double)No_C - 1.0);
  Acc_E = (E_Range[1] - E_Range[0])/((double)No_E - 1.0);
  Acc_D = 0.0; //(D_Range[1] - D_Range[0])/((double)No_D - 1.0);
  Acc_P = 0.0; //(P_Range[1] - P_Range[0])/((double)No_P - 1.0);
               //These searches here involve only colonization and extinction model 
               //parameters. The 3rd and 4th model parameters are just overwritten, 
               //in a dummy way, by the colonization and extinction values. 
  Parameter_Space_Boundaries_R_SHLIB( Space, C_Range, E_Range, C_Range, E_Range );
  Parameter_Space_Accuracies_R_SHLIB( Space, Acc_C, Acc_E, Acc_D, Acc_P );
  Parameter_Space_Initialization_R_SHLIB( Space, 
					  (*Tolerance), (*No_of_ITERATIONS), 
					  (*No_of_PARAMETERS),
                                          Index[0], Index[1], Index[2], Index[3],
                                          Discretization[0], Discretization[1], 
					  Discretization[2], Discretization[3] );
  
  Parameter_Fitting * F = (Parameter_Fitting*)calloc(1,sizeof(Parameter_Fitting));
  F->P      = P;
  F->Data   = Data;
  F->Space  = Space;
  F->Verbose = 0; // 1: Verbose // 0: Non Verbose 
  
  gsl_vector * x = gsl_vector_alloc((*No_of_PARAMETERS));
  Parameter_Model_into_Vector_Entries ( P, x,  
					Space->Parameter_Index, (*No_of_PARAMETERS) );

  if( (*Minimization) == 1 ) 
    (* Value) = GSL_Minimization_Simplex( F, x, x, 
					  &GSL_NLL_Function_Uneven );
  else if ( (*Minimization) == 0 ) 
    (* Value) = GSL_NLL_Function_Uneven( x, F);
  else {
    Rprintf(" Error in mle_Col_Ext_Uneven_Matrix_R_SHLIB(...) from\n");
    Rprintf(" file mle_Col_Ext_Uneven_Matrix_R_SHLIB.c\n");
    Rprintf(" Error in 1/0 Minimization input argument!\n ---> Minimization = %d\n", 
	   (*Minimization) ); 
  }
   Vector_Entries_into_Parameter_Model ( x, P, 
					Space->Parameter_Index, (*No_of_PARAMETERS) );
  
  (* Colonization)  = P->Colonization_Rate;     
  (* Extinction )   = P->Extinction_Rate;     
  
  P_A_R_A_M_E_T_E_R___M_O_D_E_L___F_R_E_E( P );
  
  Parameter_Space_Free(Space);
  free (Data); free (F); gsl_vector_free(x);	
}

