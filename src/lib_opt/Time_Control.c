#include <MODEL.h>

void T_I_M_E___C_O_N_T_R_O_L___A_L_L_O_C( Time_Control * Time, 
					  int OUTPUT_VARIABLES, 
					  int I_Time)
{
  int i;

  Time->Time_Vector = (double *)calloc( I_Time, sizeof(double) );

  Time->AVE = (double **)calloc(OUTPUT_VARIABLES, sizeof(double *));
  for(i = 0; i<OUTPUT_VARIABLES; i++){
    Time->AVE[i]       = (double *)calloc( I_Time, sizeof(double) );
  }

  Time->VAR = (double **)calloc(OUTPUT_VARIABLES, sizeof(double *));
  for(i = 0; i<OUTPUT_VARIABLES; i++){
    Time->VAR[i]       = (double *)calloc( I_Time, sizeof(double) );
  }

  Time->summ = (double **)calloc(OUTPUT_VARIABLES, sizeof(double *));
  for(i = 0; i<OUTPUT_VARIABLES; i++){
    Time->summ[i]       = (double *)calloc( I_Time, sizeof(double) );
  }
  
  Time->summ_var   = (double **)calloc( OUTPUT_VARIABLES, sizeof(double *) );
  for (i=0; i<OUTPUT_VARIABLES; i++){
    Time->summ_var[i]   = (double *)calloc( I_Time, sizeof(double) );
  }
}

void  T_I_M_E___C_O_N_T_R_O_L___U_P_L_O_A_D( Time_Control * Time, 
					     int I_Time, 
					     double Time_0, double Time_1 )
{
  /* Setup for the vector of sampling times */
  int i;

  Time->I_Time  = I_Time;
  Time->Time_0  = Time_0;
  Time->Time_1  = Time_1;

  Time->Delta_T = (Time->Time_1 - Time->Time_0)/(double)(I_Time-1);

  for(i=0; i<I_Time; i++){
    Time->Time_Vector[i] = Time->Time_0 + (double)i * (Time->Time_1 - Time->Time_0)/(double)(I_Time-1);
  }
  /* so that Time_Vector[0] = Time_0 and  Time_Vector[I_Time-1] = Time_1. 
     In this way, time series have I_Time points, where the first point 
     always corresponds to Time_0 and the last to Time_1;
  */
}

void T_I_M_E___C_O_N_T_R_O_L___F_R_E_E( Time_Control * Time, 
					int OUTPUT_VARIABLES )
{  
  int i;

  free (Time->Time_Vector); 

  for(i = 0; i<OUTPUT_VARIABLES; i++){
    free (Time->AVE[i]);
  }
  free (Time->AVE);

  for(i = 0; i<OUTPUT_VARIABLES; i++){
    free (Time->VAR[i]);
  }
  free (Time->VAR);
  
  for(i = 0; i<OUTPUT_VARIABLES; i++){
    free (Time->summ[i]);
  }
  free (Time->summ);

  for (i=0; i<OUTPUT_VARIABLES; i++){
    free(Time->summ_var[i]);
  }
  free (Time->summ_var);
}
