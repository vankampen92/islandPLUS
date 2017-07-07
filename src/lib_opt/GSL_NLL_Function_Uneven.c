#include <MODEL.h>

double GSL_NLL_Function_Uneven( const gsl_vector * x, void * Par )
{
  /* This GSL function allows a NLL calculation of an uneven matrix
     The input matrix is uneven, this is, each 
     presence/absence row has a different length and a 
     different, characteristic sampling vector. The matrix 
     Presence[][] has, then, the folowing structure:

      Sp
      ---
      1:  0  1  0  0  0  1  (associated to vector Sp_Time_Vector[0])   
      2:  0  1              (associated to vector Sp_Time_Vector[1]) 
      ...
      S:  1  1  0  1        (associated to vector Sp_Time_vector[S-1])
  
      The NLL is calculated and accumulated row per row.  
  */
  double NLL;
  int i,j;
  int No_of_TIMES;
  Time_Control Time;
  double ** Data; 

  Parameter_Fitting * F = (Parameter_Fitting *)Par;
  
  F->P->Time          = &Time;
  if( F->P->No_of_SPECIES != F->Data->No_of_SPECIES ) error(0,0,"Number of Species does not match: program aborted");

  int No_of_SPECIES        = F->P->No_of_SPECIES;
  double Colonization_Rate = gsl_vector_get( x, 0 );
  double Extinction_Rate   = gsl_vector_get( x, 1 );

  Data    = (double **)calloc(1, sizeof(double *) );	  
  NLL = 0.0;
  for( i=0; i<No_of_SPECIES; i++ ) {
    No_of_TIMES = F->Data->No_Sp_Time[i];

    double Time_0 = F->Data->Sp_Time[i][0];
    double Time_1 = F->Data->Sp_Time[i][No_of_TIMES-1];
    T_I_M_E___C_O_N_T_R_O_L___A_L_L_O_C( &Time, 1, No_of_TIMES );
    T_I_M_E___C_O_N_T_R_O_L___U_P_L_O_A_D( &Time, No_of_TIMES, 
					   Time_0, Time_1 );   
    Data[0] = (double *)calloc( No_of_TIMES, sizeof(double) );

    for(j=0; j<No_of_TIMES; j++) { 
      Time.Time_Vector[j] = F->Data->Sp_Time[i][j];
      Data[0][j]          = F->Data->Presence[i][j];
    }
    
    // The code for this function is located in:
    // GSL_NLLikelihood_Function.c
    NLL += NLLikelihood_Calculation( No_of_TIMES, &Time, 
				     Data, 1, 
				     Colonization_Rate, 
				     Extinction_Rate );

    T_I_M_E___C_O_N_T_R_O_L___F_R_E_E(&Time, 1);
    free(Data[0]); 
  }
  free( Data );

  return(NLL);
}

