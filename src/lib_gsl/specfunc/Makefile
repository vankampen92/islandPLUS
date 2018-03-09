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

OBJECTS =airy.o airy_der.o airy_zero.o atanint.o bessel_amp_phase.o bessel.o bessel_I0.o bessel_I1.o bessel_i.o bessel_In.o bessel_Inu.o bessel_J0.o bessel_J1.o bessel_j.o bessel_Jn.o bessel_Jnu.o bessel_K0.o bessel_K1.o bessel_k.o bessel_Kn.o bessel_Knu.o bessel_olver.o bessel_sequence.o bessel_temme.o bessel_Y0.o bessel_Y1.o bessel_y.o bessel_Yn.o bessel_Ynu.o bessel_zero.o beta.o beta_inc.o clausen.o coulomb_bound.o coulomb.o coupling.o dawson.o debye.o dilog.o elementary.o ellint.o elljac.o erfc.o exp.o expint3.o expint.o fermi_dirac.o gamma.o gamma_inc.o gegenbauer.o hermite.o hyperg_0F1.o hyperg_1F1.o hyperg_2F0.o hyperg_2F1.o hyperg.o hyperg_U.o laguerre.o lambert.o   log.o mathieu_angfunc.o mathieu_charv.o mathieu_coeff.o mathieu_radfunc.o mathieu_workspace.o poch.o pow_int.o psi.o result.o shint.o sinint.o synchrotron.o transport.o trig.o zeta.o

.c.o:
	${CC} $(CFLAGS) ${EXTRA_INCLUDE_PATH} -c $<

all: ${OBJECTS} 

clean:  
	rm -f *.o *~
