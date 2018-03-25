EXTRA_INCLUDE_PATH = -I./lib_opt/ -I../include/ -I../lib/ -I/usr/share/R/include/ -I/Library/Frameworks/R.framework/Resources/include/
## R_SHLIB is my own compilation FLAG for the generation of libraries compatible 
## with the generation of a R-callable shared objects.  
RSH =R_SHLIB

CFLAGS+=-DHAVE_INLINE # -D${RSH} ${EXTRA_INCLUDE_PATH} 
###CFLAGS= -Wall -O21

CC=gcc -fPIC -DSHARED_OBJECT 

RANLIB=ranlib 

OBJECTS=Time_Control.o Parameter_Model.o Parameter_Space_R_SHLIB.o SP_Matrix_Data.o Parameter_Index_Checking_Order.o Checking_for_Parameter_Correctness.o GSL_NLLikelihood_Function.o GSL_NLL_Function_Uneven.o GSL_MacKenzie_NLLikelihood_Function.o GSL_MacKenzie_NLL_Uneven_Function.o MacKenzie_NLLikelihood_Function.o GSL_Minimization_Simplex.o mle_Col_Ext_Uneven_Matrix_R_SHLIB.o MODEL_SELECTION_UPGMA_R_SHLIB.o mle_MacKenzie_Uneven_Matrix_R_SHLIB.o MODEL_SELECTION_UPGMA_MacKENZIE_R_SHLIB.o

.SUFIXES: .o .c

.c.o: 
	${CC} $(CFLAGS) -DHAVE_INLINE -D${RSH} ${EXTRA_INCLUDE_PATH} -c $<

IBT_Optimization: ${OBJECTS} ../include/MODEL.h
	${AR} -q libda_IBT_Optimization.a ${OBJECTS}
	${RANLIB} libda_IBT_Optimization.a

clean:  
	rm -f *.o *.a *~

