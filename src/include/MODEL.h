#define MODEL_BASIC_LIBRARIES
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

// #include <complex.h>

#include <time.h>
#include <assert.h>
#include <gsl/gsl_sf.h> 
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_multimin.h> 
#include <gsl/gsl_histogram.h>

// Curiously R.h is incompatible with the C general header complex.h
// If you are not requiring funcions on the complex plane, you are OK!!! 
#include <R.h>
 
#include "MODEL_Time_Control_STRUCT_DEF.h"
#include "MODEL_Parameter_Model_STRUCT_DEF.h"
#include "MODEL_SP_Matrix_Data_STRUCT_DEF.h"
#include "MODEL_Parameter_Space_STRUCT_DEF.h"
#include "MODEL_Parameter_Fitting_STRUCT_DEF.h"

#include "../lib/da_IBT_Functions.h"
#include "../lib_opt/da_IBT_Optimization.h"

#define MAX(A,B) ((A)>(B)? A:B)
#define MIN(A,B) ((A)<(B)? A:B)

#define No_of_SPECIES_MAX 10000
#define MODEL_PARAMETERS_MAXIMUM 4
#define No_of_TIME_OBS_MAX 100000 // Maximum number of temporal observations 
                                  // (total temporal columns in the data matrix) 
#define ROW_TAG_LENGTH_MAX 100    // Maximum number of characters of the label
                                  // of each row of a given Presence Absence 
                                  // data

