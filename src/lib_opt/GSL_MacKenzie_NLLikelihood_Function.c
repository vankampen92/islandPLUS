#include <MODEL.h>

double GSL_MacKenzie_NLLikelihood_Function ( const gsl_vector * x, void * Par )
{
  Parameter_Fitting * F     = (Parameter_Fitting *)Par;
  Parameter_Space   * Space = F->Space; 
  Parameter_Model   * P     = F->P; 

  if ( F->P->No_of_SPECIES != F->Data->No_of_SPECIES )
    error(0,0, "Number of Species does not match: program aborted");

  int * Transects     = F->Data->Transects;
  double * T          = F->Data->Time_Vector;
  int No_of_SPECIES   = F->P->No_of_SPECIES;
  int n               = F->Data->No_of_TIMES;
  int N               = F->Data->Total_No_of_TRANSECTS;

  Vector_Entries_into_Parameter_Model ( x, P, 
					Space->Parameter_Index,
					Space->No_of_PARAMETERS );

  double Colonization_Rate   = P->Colonization_Rate;
  double Extinction_Rate     = P->Extinction_Rate;
  double Detectability_Value = P->Detectability_Value;
  double Phi_0               = P->Phi_0;

  double ** Data = F->Data->Presence;
  
  if (F->P->RATES == 0) { 
    // Input parameters are transition probabilities 
    // and should be converted into true rates before 
    // calling the neg likelihood evaluation function
    Probability_Rates( Colonization_Rate, Extinction_Rate,
		       &Colonization_Rate, &Extinction_Rate,
		       1.0 ) ;
    // Intersampling Time is assumed to be 1.0.
  }

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

  double NLL = MacKenzie_NLLikelihood_Calculation( Data, No_of_SPECIES, N,
						   T, Transects, n, 
						   Colonization_Rate, 
						   Extinction_Rate,
						   Detectability_Value,
						   Phi_0 );
  return(NLL);
}
