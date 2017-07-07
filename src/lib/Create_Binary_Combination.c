/*
   This code creates an ordered list of binary numbers:
   from 1 to 2^{MAX_LENGTH}. Binary numbers are stored 
   in a int ** Binary_Combination array:
	
	0:                  Binary_Combination[0] = {00000000}
	1:                  Binary_Combination[1] = {00000001}
	2:                  Binary_Combination[2] = {00000010}
	.
.	.
	.
	2^{MAX_LENGTH}-1    Binary_Combination[2^{MAX_LENGTH}-1]
   Compilation: 
   $ gcc -o BINARY_LIST main.c -lm

   Execution: 
   $ ./BINARY_LIST 10
*/
#include "HEADERS.h"
#include "da_IBT_Functions.h"

void int_buffer_rec(int ** Number_List, int N, 
		    int * number, int n, int length) 
{
    int i;	
    static int m = 0;

    if(n > 0) {
        number[length - n] = 0;
        int_buffer_rec(Number_List, N, number, n - 1, length);
        number[length - n] = 1;
        int_buffer_rec(Number_List, N, number, n - 1, length);
    }
    else { 
        for(i = 0; i < length; ++i) {
            // Rprintf("%u", number[i]);
	    Number_List[m][i]=number[i];
        }
        // Rprintf("\n");
        m++;
    }

    if ( m == N ) m = 0;
}

void Create_Binary_Combination( int ** Binary_Combination, int N, 
				int LENGTH )
{
  /* This function creates exhaustively all 1/0 strings of a given 
     length LENGTH.      

     The numbers of strings to create should be given as the N input
     parameter. Enough space should have been reserved (in the parent function) in 
     Binary_Combination[][] to store all of them */

                int * number  = (int *)calloc(LENGTH, sizeof( int ) );

                int_buffer_rec(Binary_Combination, N, number, LENGTH, LENGTH);

                free(number);
}
