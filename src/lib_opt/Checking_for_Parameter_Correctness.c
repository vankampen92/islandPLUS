#include <MODEL.h>

int Checking_for_Parameter_Correctness( double Colonization_Rate, 
					double Extinction_Rate,
					double Detectability_Value, 
					double Phi_0 )
{
  /* This little function is required by the procedures evaluating a
     likelihood. These procedures checked if the parameters have taken
     non-sense values, such as negative values, etc, at any minimization
     iteration */

  /* Model parameters are correctly bounded if PAR_DEFINITION 
     flag is left unchanged and returned as 1. 
  */
  int PAR_DEFINITION;
  
  PAR_DEFINITION = 1;  

  if(Colonization_Rate <= 0.0) PAR_DEFINITION = 0;

  if(Extinction_Rate <= 0.0) PAR_DEFINITION = 0;

  if(Detectability_Value <= 0 || Detectability_Value >= 1 ) PAR_DEFINITION = 0;

  if(Phi_0 <= 0 || Phi_0 >= 1) PAR_DEFINITION = 0;
  
  return (PAR_DEFINITION);
}
