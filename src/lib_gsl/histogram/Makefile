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

OBJECTS =  add2d.o add.o calloc_range2d.o calloc_range.o copy2d.o copy.o file2d.o file.o get2d.o get.o init2d.o init.o maxval2d.o maxval.o oper2d.o oper.o params2d.o params.o pdf2d.o pdf.o reset2d.o reset.o stat2d.o stat.o urand.o 

.c.o:
	${CC} $(CFLAGS) ${EXTRA_INCLUDE_PATH} -c $<

all: 	${OBJECTS}

clean:  
	rm -f *.o *~
