EXTRA_INCLUDE_PATH = -I./lib/ -I../include/ -I/usr/share/R/include/

CFLAGS+=-DHAVE_INLINE # ${EXTRA_INCLUDE_PATH} 

CC=gcc -fPIC -DSHARED_OBJECT

RANLIB=ranlib

OBJECTS=UPGMA_CLUSTERING_PARTITION.o upgma_clustering.o EUCLIDIAN_DISTANCE_MATRIX.o Probability_Rates.o Transition_Matrix.o Rates_into_Probabilities.o IO_Missing_Values.o IO_AKAIKE_Model_Selection.o Counting_Type_of_Transitions.o Create_Binary_Combination.o


.SUFIXES:       .o .c

.c.o: 
	${CC} $(CFLAGS) -DHAVE_INLINE ${EXTRA_INCLUDE_PATH} -c $<

IBT_Functions:	${OBJECTS} HEADERS.h
	${AR} -q libda_IBT_Functions.a ${OBJECTS}
	${RANLIB} libda_IBT_Functions.a

clean:  
	rm -f *.o *.a *~


