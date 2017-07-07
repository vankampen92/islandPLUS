#include <MODEL.h>

double GSL_MacKenzie_NLL_Uneven_Function ( const gsl_vector * x, void * Par )
{
  /* This function evaluates MacKenzie likelihood on a per row basis. 
     This is so required because the time and transect structure of the data matrix 
     differs from row to row. When data matrix keeps the same time structure for all 
     rows, the equivalent likelihood function to call is: 
 
     GSL_MacKenzie_NLLikelihood_Function(const gsl_vector * x, void * Par); 
  */ 
  int i,j;
  double NLL;
  int S,n,N;
  int * Transects;
  double * Time;
  
  Parameter_Fitting * F     = (Parameter_Fitting *)Par;
  Parameter_Space   * Space = F->Space; 
  Parameter_Model   * P     = F->P;
  
  int No_of_SPECIES   = F->P->No_of_SPECIES;

  if ( F->P->No_of_SPECIES != F->Data->No_of_SPECIES ) 
     error(0,0, "Number of Species does not match: program aborted");
  
  Vector_Entries_into_Parameter_Model ( x, P, 
					Space->Parameter_Index, 
					Space->No_of_PARAMETERS );

  double Colonization_Rate   = P->Colonization_Rate;
  double Extinction_Rate     = P->Extinction_Rate;
  double Detectability_Value = P->Detectability_Value;
  double Phi_0               = P->Phi_0;

  if (F->P->RATES == 0) { 
      // Input parameters are transition probabilities 
      // and should be converted into true rates before 
      // calling the neg likelihood evaluation function
      Probability_Rates( Colonization_Rate, Extinction_Rate,
			 &Colonization_Rate, &Extinction_Rate,
			 1.0 ) ;
      // Intersampling Time is assumed to be 1.0.
  }
  
  NLL = 0.0;
  for( i= 0; i<No_of_SPECIES; i++ ) { 

    Transects  = F->Data->Sp_Transects[i];
    Time       = F->Data->Sp_Time[i];
    S          = 1; 
    n          = F->Data->No_Sp_Time[i];
    N          = F->Data->Sp_Total_No_Transects[i];

    double ** Data = (double **)calloc(1, sizeof(double *) );
    Data[0]        = (double *)calloc(N, sizeof(double) );
    for( j = 0; j<N; j++ ) Data[0][j] = F->Data->Presence[i][j];

    /* The code corresponding to:
     *
     *              MacKenzie_NLLikelihood_Calculation(...)
     *
     * is in the file: 
     *  
     *              MacKenzie_NLLikelihood_Function.c 
     *
     * where the whole likelikehood evaluation is actually performed!!!  
     */
    
    NLL += MacKenzie_NLLikelihood_Calculation( Data, S, N,
					       Time, Transects, n, 
					       Colonization_Rate, 
					       Extinction_Rate,
					       Detectability_Value,
					       Phi_0 );
    free(Data[0]); 
    free(Data); 
  }
  
  return(NLL);
}
