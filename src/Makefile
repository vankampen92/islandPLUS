###################################################################
#       B E G I N :  Build shared library of functions to be callable from R ('isandPLUS.so')
SOURCES=MODEL_SELECTION_UPGMA_R_FUNCTION.c MODEL_SELECTION_UPGMA_MacKENZIE_R_FUNCTION.c R_SHLIB___mle_MacKenzie_NLL_Uneven_Minimization.c R_SHLIB___mle_NLLikelihood_Minimization.c R_SHLIB___mle_MacKenzie_NLLikelihood_Minimization.c Counting_Replicates_per_Time.c
SHLIB=islandPLUS.so
AR=libtool -static -arch_only x86_64 
#
all:  	IBT_Functions IBT_Optimization IBT_gsl
	echo "Building shared object "${SHLIB}
	R CMD SHLIB -o ${SHLIB} ${SOURCES} 
# Auxiliary libraries 
IBT_gsl: 
	echo "Compiling Library da_gsl.a needed by shared object "${SHLIB}
	cd ./lib_gsl; ${MAKE} MAKE=${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)" AR="$(AR)" RANLIB="$(RANLIB)"
IBT_Functions: 
	echo "Compiling Library da_IBT_Functions.a needed by shared object "${SHLIB}
	cd ./lib; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)" AR="$(AR)" RANLIB="$(RANLIB)"
IBT_Optimization: 
	echo "Compiling Library da_IBT_Optimization.a needed by shared object "${SHLIB}
	cd ./lib_opt; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)" AR="$(AR)" RANLIB="$(RANLIB)"
##################################################################
clean: 
	rm *.o; cd ./gsl; rm -r gsl-2.4; cd ./lib; ${MAKE} clean; cd ../lib_opt; ${MAKE} clean; cd ../lib_gsl; ${MAKE} clean
#       E N D :  Example to build an R shared library from C code
###################################################################
