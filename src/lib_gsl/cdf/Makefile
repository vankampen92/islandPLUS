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

OBJECTS =beta.o beta_inc.o betainv.o binomial.o cauchy.o cauchyinv.o chisq.o chisqinv.o exponential.o exponentialinv.o exppow.o fdist.o fdistinv.o flat.o flatinv.o gamma.o gammainv.o gauss.o gaussinv.o geometric.o gumbel1.o gumbel1inv.o gumbel2.o gumbel2inv.o hypergeometric.o laplace.o laplaceinv.o logistic.o logisticinv.o lognormal.o lognormalinv.o nbinomial.o pareto.o paretoinv.o pascal.o poisson.o rayleigh.o rayleighinv.o tdist.o tdistinv.o weibull.o weibullinv.o 

.c.o:
	${CC} $(CFLAGS) ${EXTRA_INCLUDE_PATH} -c $<

all: ${OBJECTS} 

clean:  
	rm -f *.o *~
