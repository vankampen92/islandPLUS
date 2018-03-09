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

OBJECTS =bernoulli.o beta.o bigauss.o binomial.o binomial_tpe.o cauchy.o chisq.o dirichlet.o discrete.o erlang.o exponential.o exppow.o fdist.o flat.o gamma.o gauss.o gausstail.o gausszig.o geometric.o gumbel.o hyperg.o landau.o laplace.o levy.o logarithmic.o logistic.o lognormal.o multinomial.o mvgauss.o nbinomial.o pareto.o pascal.o poisson.o rayleigh.o shuffle.o sphere.o tdist.o weibull.o

.c.o:
	${CC} $(CFLAGS) ${EXTRA_INCLUDE_PATH} -c $<

all: ${OBJECTS} 

clean:  
	rm -f *.o *~
