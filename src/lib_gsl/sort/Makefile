# ==========================================
# CUSTUMIZED GSL-based UTILITY LIBRARY make
# ==========================================
#
# ==========================================
# Gnu C make to build a gsl limited Library 
# The generated library is libda_gsl.a 
# ==========================================
EXTRA_INCLUDE_PATH = -I./ -I../../include/ -I../../include/gsl/ -I/usr/share/R/include/
CFLAGS=-g -DHAVE_INLINE 
###CFLAGS= -Wall -O21 
CC=gcc -fPIC -DSHARED_OBJECT

OBJECTS = sort.o sortind.o sortvec.o sortvecind.o subset.o subsetind.o  

.c.o:
	${CC} $(CFLAGS) ${EXTRA_INCLUDE_PATH} -c $<

all: 	${OBJECTS}

clean:  
	rm -f *.o *~
