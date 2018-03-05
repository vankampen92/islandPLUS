EXTRA_INCLUDE_PATH = -I./lib_gsl/ -I../include/ -I../lib/ -I/usr/share/R/include/
## R_SHLIB is my own compilation FLAG for the generation of libraries compatible 
## with the generation of a R-callable shared objects.  
RSH =R_SHLIB

CFLAGS+=-DHAVE_INLINE # -D${RSH} ${EXTRA_INCLUDE_PATH} 
###CFLAGS= -Wall -O21

CC=gcc -fPIC -DSHARED_OBJECT 

RANLIB=ranlib

OBJECTS=

.SUFIXES: .o .c

.c.o: 
	${CC} $(CFLAGS) -DHAVE_INLINE -D${RSH} ${EXTRA_INCLUDE_PATH} -c $<

IBT_Optimization: ${OBJECTS} ../include/MODEL.h
	${AR} -q libda_IBT_Optimization.a ${OBJECTS}
	${RANLIB} libda_IBT_Optimization.a

clean:  
	rm -f *.o *.a *~

