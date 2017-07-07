#include <MODEL.h>

void Counting_Replicates_per_Time(double * Vector, int * Temporal_Observations,
				  double * Time_Vector, int * Transects, 
				  int * No_of_TIMES)
{
  /* This functions counts transects/replicates per sampling time. The number 
     of replicates per sampling time can perfectly differ. 
     
     Input:
     . Vector[], array of total sampling times (including time repetitions);
     . Temporal_Observations, total number of observations (length of Vector[])
     
     Output:
     . Time_Vector, array of single sampling times (with no time repetitions)
     . Transects, array counting the number of repetitions or replicates per 
     sampling time. 
     . No_of_TIMES, true number of different times (length of both Transects and 
     Time_Vector arrays)
  */
  double y; 
  int nt, m, n;
  int N; 
  
  N = (*Temporal_Observations);
  nt = 0;
  m  = 0;
  n  = 0;
  y=Vector[0];
  while (n < N) {
   if( Vector[n] == y ){
     Time_Vector[m] = y;
     nt++;
     n++;
   }
   else {
     y = Vector[n];
     Transects[m++] = nt;
     nt = 0;
   }
  }
  Transects[m] = nt;
  (* No_of_TIMES) = 1 + m;
}
