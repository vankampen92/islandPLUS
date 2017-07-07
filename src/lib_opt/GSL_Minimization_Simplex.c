#include <MODEL.h>

/* #define TOLERANCE 1.0E-08        */
/* #define MAX_No_of_ITERATIONS 100 */

double GSL_Minimization_Simplex (Parameter_Fitting * F, 
				 gsl_vector * Initial_Guess, 
				 gsl_vector * Solution, 
				 double ( * Function )( const gsl_vector * , void * ) )
{
  int i;
  //int key;
  double value;

  Parameter_Space * Space  = F->Space;
  // Parameter_Model * P      = F->P;

  const gsl_multimin_fminimizer_type *T = 
    gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss; 
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  int No_of_PARAMETERS     = Space->No_of_PARAMETERS;
  double TOLERANCE         = Space->TOLERANCE; 
  int MAX_No_of_ITERATIONS = Space->MAX_No_of_ITERATIONS; 
  
  /* Set initial step sizes */
  ss = gsl_vector_alloc ( No_of_PARAMETERS );
  // gsl_vector_set_all (ss, 0.01);
  gsl_vector_memcpy (ss, Space->Accuracy );

  /* Initialize method and iterate */
 
  if ( F->Verbose == 1 ) Rprintf("No_of_PARAMETERS = %d\n", No_of_PARAMETERS); 

  minex_func.n = No_of_PARAMETERS;
  minex_func.f = Function;
  minex_func.params = (void *)F;

  //s = gsl_multimin_fminimizer_alloc (T, Space->No_of_PARAMETERS);
  s = gsl_multimin_fminimizer_alloc (T, No_of_PARAMETERS);
  gsl_multimin_fminimizer_set (s, &minex_func, Initial_Guess, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if (status) 
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, TOLERANCE);

       if (F->Verbose == 1) { 
	 if (status == GSL_SUCCESS)
	   {
	     Rprintf("converged to minimum at\n");
	   }
     
	 Rprintf("%5d ", (int)iter);
	 // for( i=0; i<Space->No_of_PARAMETERS; i++ ) {
	 // key  = Space->Parameter_Index[i];
	 // Rprintf("%s = %10.3e\t", P->Symbol[key], gsl_vector_get (s->x, i) ); 
	 for( i=0; i<No_of_PARAMETERS; i++ ) {
	   if (i == 0) Rprintf("Colonization = %10.3e; ", gsl_vector_get (s->x, i) ); 
	   if (i == 1) Rprintf("Extinction   = %10.3e; ", gsl_vector_get (s->x, i) ); 
	   if (i == 2) Rprintf("Detectability = %10.3e; ", gsl_vector_get (s->x, i) ); 
	   if (i == 3) Rprintf("Phi_0   = %10.3e\t", gsl_vector_get (s->x, i) ); 
	 }
	 Rprintf("f() = %7.3f size = %.3f\n", s->fval, size);
       }

    }
  while (status == GSL_CONTINUE && iter < MAX_No_of_ITERATIONS );
     
  gsl_vector_memcpy (Solution, s->x ); 
  value = s->fval; // Min Value !!! 

  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  //return status;
  return( value );
}
