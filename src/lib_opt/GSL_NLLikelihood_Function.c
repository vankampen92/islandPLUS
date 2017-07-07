#include <MODEL.h>

#define CHECKING_MATRIX_ENTRIES
// When the CHECKING_MATRIX_ENTRIES flag is defined 
// the presence data matrix has to have only 0/1 entries.  
// Otherwise, the likelihood calculation fails.
// If there are missing data, and some entries are 
// not 0/1, the likelihood calculation can still proceed
// if this flag is not defined (the "define" line above is 
// commented out)  
 
double GSL_NLLikelihood_Function ( const gsl_vector * x, void * Par )
{
  Parameter_Fitting * F = (Parameter_Fitting *)Par;

  if( F->P->No_of_SPECIES != F->Data->No_of_SPECIES ) error(0,0,"Number of Species does not match: program aborted");

  Time_Control * T    = F->P->Time;
  int No_of_SPECIES   = F->P->No_of_SPECIES;
  int n               = F->P->Time->I_Time;

  double Colonization_Rate = gsl_vector_get( x, 0 );
  double Extinction_Rate   = gsl_vector_get( x, 1 );

  double ** Data = F->Data->Presence;
  
  double NLL = NLLikelihood_Calculation( n, T, 
					 Data, 
					 No_of_SPECIES,
					 Colonization_Rate, 
					 Extinction_Rate );
  return(NLL);
}

double NLLikelihood_Calculation ( int n, Time_Control * Time, 
				  double ** Presence_Data, 
				  int No_of_SPECIES,
				  double Colonization_Rate, 
				  double Extinction_Rate )
{
  /* Input arguments:
     . n : Number of empirical observations (sampling times)
     . Time : Structure containing sampling times. 
     . Presence_Data[][] : 1/0 Data matrix. This matrix may have some 
      missing values, which are always maked with the label, the double
      number '0.1'. If some other level is used. This code will fail. 
     . No_of_SPECIES: Number of rows of Data matrix. 
     . Colonization_Rate: colonization rate in T^{-1} units
     . Extintion_Rate: extinction rate in T^{-1} units

     Output arguments:
     . Neg Log Likelihood 
     (see App A: "Likelihood Estimation of Transition Probabilities", Eq A8)
  */ 
  int i,j;
  int k_0, k_1;
  double P_k_0, P_k_1;
  double X, XX, NLL;
  
  int PAR_DEFINITION = Checking_for_Parameter_Correctness( Colonization_Rate, 
							   Extinction_Rate,
							   0.5, 
							   0.5 );	
					
  if ( PAR_DEFINITION  ==  0 ) {
   
    X = DBL_MAX;
   
    return (X);
  } 
  else {

    double ** TM = (double **)calloc(2, sizeof(double *) );
    TM[0]= (double *)calloc(2, sizeof(double) );
    TM[1]= (double *)calloc(2, sizeof(double) );
    
    X = 0.0; NLL = 0.0;
    for(i=1; i<n; i++) {
      
      Time->Delta_T = Time->Time_Vector[i] - Time->Time_Vector[i-1];
      Transition_Matrix( TM, 2,2,
			 Colonization_Rate,
			 Extinction_Rate,
			 Time->Delta_T );
      
      for(j=0; j<No_of_SPECIES; j++) {
	/* Transiction: k_0 ---> k_1 */
	P_k_0 = Presence_Data[j][i-1];
	P_k_1 = Presence_Data[j][i];
	
	if ( (P_k_0==0.0 || P_k_0==1.0) && (P_k_1==0.0 || P_k_1==1.0) ) {
	  
	  k_0 = (int)Presence_Data[j][i-1];
	  k_1 = (int)Presence_Data[j][i];
	  
          if( TM[k_1][k_0] > 0.0 && TM[k_1][k_0] <= 1.0) XX = log(TM[k_1][k_0]);
	  else { NLL = DBL_MAX; return(NLL); }

	  X = -XX;
	  
	}
	else {
#if defined CHECKING_MATRIX_ENTRIES
	  Rprintf(" Error in initial Presence Data\n");
	  Rprintf(" when evaluation Neg LogLikelihood\n");
	  Rprintf(" in function GSL_NLLikelihood_Function(...)\n");
	  Rprintf(" Some matrix entries are not either 0 or 1\n");
	  Rprintf(" (see GSL_NLLikelihood_Function.c)\n");
	  Rprintf(" The program will exit\n");
	  error(0, 0, "Program has aborted");
#endif
	}
        NLL += X;	
      }   
    }
    
    free(TM[0]); free(TM[1]);
    free(TM);

    return(NLL);
  }
}

double NLLikelihood_Calculation_Transition_Probabilities ( int n, Time_Control * Time, 
							   double ** Presence_Data, 
							   int No_of_SPECIES,
							   double Colonization_Probability, 
							   double Extinction_Probability )
{
  /* Input arguments:
     . n : Number of empirical observations
     . Time : Structure containing sampling times. 
     . Presence_Data[][] : 1/0 Data matrix. This matrix may have some 
      missing values, which are always maked with the label, the double
      number '0.1'. If some other level is used. This code will fail. 
     . No_of_SPECIES: Number of rows of Data matrix. 
     . Colonization: colonization transition probability calculated
     with the first n_Time columns (n_Time < n);
     . Extintion: Extinction transition probability calculated
     with the first n_Time columns (n_Time < n);
     
     Output arguments:
     . Neg Log Likelihood 
     (see App A: "Likelihood Estimation of Transition Probabilities", Eq A8)
  */
  int i,j;
  int k_0, k_1;
  double P_k_0, P_k_1;
  double X;
  double Colonization_Rate, Extinction_Rate;

  double ** TM = (double **)calloc(2, sizeof(double *) );
  TM[0]= (double *)calloc(2, sizeof(double) );
  TM[1]= (double *)calloc(2, sizeof(double) );

  Probability_Rates( Colonization_Probability, Extinction_Probability,
                     &Colonization_Rate, &Extinction_Rate,
                     1.0 ) ;
  X = 0.0;
  for(i=1; i<n; i++) {

    Time->Delta_T = Time->Time_Vector[i] - Time->Time_Vector[i-1];
    Transition_Matrix( TM, 2,2,
                       Colonization_Rate,
                       Extinction_Rate,
                       Time->Delta_T );
    
    for(j=0; j<No_of_SPECIES; j++) {
      /* Transiction: k_0 ---> k_1 */
      P_k_0 = Presence_Data[j][i-1];
      P_k_1 = Presence_Data[j][i];
      
      if  (P_k_0 == 0.0 || P_k_0 == 1.0 ) {
	if  (P_k_1 == 0.0 || P_k_1 == 1.0 ) {
	  
	  k_0 = (int)Presence_Data[j][i-1];
	  k_1 = (int)Presence_Data[j][i];
	  
	  X += log(TM[k_1][k_0]);
	  
	}
	else {
	  if( P_k_1 != 0.1 ) error(0,0,"Missing Value Problem: program aborted");
	}	  
      }
      else {
	  if( P_k_0 != 0.1 ) error(0,0,"Missing Value Problem: program aborted");
      }
    }
  }

  free(TM[0]); free(TM[1]);
  free(TM);
  
  X = -X;
  return(X);
}



