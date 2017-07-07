typedef struct Time_Controlinfo
{
#include "include.Time_Control.global.h" /* Control variable I_Time */

  double * Time_Vector;

  double ** summ;
  double ** summ_var;

  double ** AVE;
  double ** VAR;
  
}Time_Control;

// These functions are compiled in the library:
// libda_IBT_Optimization.a
// (see ./lib_opt/Time_Control.c)

void T_I_M_E___C_O_N_T_R_O_L___A_L_L_O_C( Time_Control * , int , int );

void  T_I_M_E___C_O_N_T_R_O_L___U_P_L_O_A_D( Time_Control * , 
					     int, 
					     double, double );

void T_I_M_E___C_O_N_T_R_O_L___F_R_E_E( Time_Control * Time, int);

