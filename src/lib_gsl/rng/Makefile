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

OBJECTS =borosh13.o cmrg.o coveyou.o default.o file.o fishman18.o fishman20.o fishman2x.o gfsr4.o inline.o knuthran2002.o knuthran2.o knuthran.o lecuyer21.o minstd.o mrg.o mt.o r250.o ran0.o ran1.o ran2.o ran3.o rand48.o rand.o random.o randu.o ranf.o ranlux.o ranlxd.o ranlxs.o ranmar.o rng.o schrage.o slatec.o taus113.o taus.o transputer.o tt.o types.o uni32.o uni.o vax.o waterman14.o zuf.o 

.c.o:
	${CC} $(CFLAGS) ${EXTRA_INCLUDE_PATH} -c $<

all: ${OBJECTS}

clean:  
	rm -f *.o *~
