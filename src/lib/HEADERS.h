#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <error.h>

#define MAX(A,B) ((A)>(B)? A:B)
#define MIN(A,B) ((A)<(B)? A:B)

#include <time.h>
#include <assert.h>
#include <limits.h>

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

#include <R.h>

