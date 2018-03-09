# ==========================================
# CUSTUMIZED GSL-based UTILITY LIBRARY make
# ==========================================
#
# ==========================================
# Gnu C make to build a gsl limited Library 
# The generated library is libda_gsl.a 
# ==========================================
MAKE=make
CFLAGS=-g -DHAVE_INLINE
CC=gcc -fPIC -DSHARED_OBJECT

GSL_LIBRARY=./
GSL_BLOCK=${GSL_LIBRARY}block/
GSL_VECTOR=${GSL_LIBRARY}vector/
GSL_ERR=${GSL_LIBRARY}err/
GSL_RNG=${GSL_LIBRARY}rng/
GSL_RANDIST=${GSL_LIBRARY}randist/
GSL_CDF=${GSL_LIBRARY}cdf/	
GSL_SPECFUNC=${GSL_LIBRARY}specfunc/
GSL_MATRIX=${GSL_LIBRARY}matrix/
GSL_SORT=${GSL_LIBRARY}sort/
GSL_PERMUTATION=${GSL_LIBRARY}permutation/
GSL_MULTIMIN=${GSL_LIBRARY}multimin/
GSL_HISTOGRAM=${GSL_LIBRARY}histogram/
GSL_SYS=${GSL_LIBRARY}sys/
GSL_BLAS=${GSL_LIBRARY}blas/
RANLIB=ranlib

OBJECTS = ./vector/copy.o ./vector/file.o ./vector/init.o ./vector/minmax.o ./vector/oper.o ./vector/prop.o ./vector/reim.o ./vector/subvector.o ./vector/swap.o ./vector/vector.o ./vector/view.o ./matrix/copy.o ./matrix/file.o ./matrix/getset.o ./matrix/init.o ./matrix/matrix.o ./matrix/minmax.o  ./matrix/oper.o  ./matrix/prop.o  ./matrix/rowcol.o ./matrix/submatrix.o  ./matrix/swap.o ./matrix/view.o ./err/error.o ./err/message.o ./err/stream.o ./err/strerror.o ./rng/borosh13.o ./rng/cmrg.o ./rng/coveyou.o ./rng/default.o ./rng/file.o ./rng/fishman18.o ./rng/fishman20.o ./rng/fishman2x.o ./rng/gfsr4.o ./rng/inline.o ./rng/knuthran2002.o ./rng/knuthran2.o ./rng/knuthran.o ./rng/lecuyer21.o ./rng/minstd.o ./rng/mrg.o ./rng/mt.o ./rng/r250.o ./rng/ran0.o ./rng/ran1.o ./rng/ran2.o ./rng/ran3.o ./rng/rand48.o ./rng/rand.o ./rng/random.o ./rng/randu.o ./rng/ranf.o ./rng/ranlux.o ./rng/ranlxd.o ./rng/ranlxs.o ./rng/ranmar.o ./rng/rng.o ./rng/schrage.o ./rng/slatec.o ./rng/taus113.o ./rng/taus.o ./rng/transputer.o ./rng/tt.o ./rng/types.o ./rng/uni32.o ./rng/uni.o ./rng/vax.o ./rng/waterman14.o ./rng/zuf.o ./cdf/beta.o ./cdf/beta_inc.o ./cdf/betainv.o ./cdf/binomial.o ./cdf/cauchy.o ./cdf/cauchyinv.o ./cdf/chisq.o ./cdf/chisqinv.o ./cdf/exponential.o ./cdf/exponentialinv.o ./cdf/exppow.o ./cdf/fdist.o ./cdf/fdistinv.o ./cdf/flat.o ./cdf/flatinv.o ./cdf/gamma.o ./cdf/gammainv.o ./cdf/gauss.o ./cdf/gaussinv.o ./cdf/geometric.o ./cdf/gumbel1.o ./cdf/gumbel1inv.o ./cdf/gumbel2.o ./cdf/gumbel2inv.o ./cdf/hypergeometric.o ./cdf/laplace.o ./cdf/laplaceinv.o ./cdf/logistic.o ./cdf/logisticinv.o ./cdf/lognormal.o ./cdf/lognormalinv.o ./cdf/nbinomial.o ./cdf/pareto.o ./cdf/paretoinv.o ./cdf/pascal.o ./cdf/poisson.o ./cdf/rayleigh.o ./cdf/rayleighinv.o ./cdf/tdist.o ./cdf/tdistinv.o ./cdf/weibull.o ./cdf/weibullinv.o ./specfunc/airy.o ./specfunc/airy_der.o ./specfunc/airy_zero.o ./specfunc/atanint.o ./specfunc/bessel_amp_phase.o ./specfunc/bessel.o ./specfunc/bessel_I0.o ./specfunc/bessel_I1.o ./specfunc/bessel_i.o ./specfunc/bessel_In.o ./specfunc/bessel_Inu.o ./specfunc/bessel_J0.o ./specfunc/bessel_J1.o ./specfunc/bessel_j.o ./specfunc/bessel_Jn.o ./specfunc/bessel_Jnu.o ./specfunc/bessel_K0.o ./specfunc/bessel_K1.o ./specfunc/bessel_k.o ./specfunc/bessel_Kn.o ./specfunc/bessel_Knu.o ./specfunc/bessel_olver.o ./specfunc/bessel_sequence.o ./specfunc/bessel_temme.o ./specfunc/bessel_Y0.o ./specfunc/bessel_Y1.o ./specfunc/bessel_y.o ./specfunc/bessel_Yn.o ./specfunc/bessel_Ynu.o ./specfunc/bessel_zero.o ./specfunc/beta.o ./specfunc/beta_inc.o ./specfunc/clausen.o ./specfunc/coulomb_bound.o ./specfunc/coulomb.o ./specfunc/coupling.o ./specfunc/dawson.o ./specfunc/debye.o ./specfunc/dilog.o ./specfunc/elementary.o ./specfunc/ellint.o ./specfunc/elljac.o ./specfunc/erfc.o ./specfunc/exp.o ./specfunc/expint3.o ./specfunc/expint.o ./specfunc/fermi_dirac.o ./specfunc/gamma.o ./specfunc/gamma_inc.o ./specfunc/gegenbauer.o ./specfunc/hermite.o ./specfunc/hyperg_0F1.o ./specfunc/hyperg_1F1.o ./specfunc/hyperg_2F0.o ./specfunc/hyperg_2F1.o ./specfunc/hyperg.o ./specfunc/hyperg_U.o ./specfunc/laguerre.o ./specfunc/lambert.o   ./specfunc/log.o ./specfunc/mathieu_angfunc.o ./specfunc/mathieu_charv.o ./specfunc/mathieu_coeff.o ./specfunc/mathieu_radfunc.o ./specfunc/mathieu_workspace.o ./specfunc/poch.o ./specfunc/pow_int.o ./specfunc/psi.o ./specfunc/result.o ./specfunc/shint.o ./specfunc/sinint.o ./specfunc/synchrotron.o ./specfunc/transport.o ./specfunc/trig.o ./specfunc/zeta.o ./randist/bernoulli.c ./randist/beta.c ./randist/bigauss.c ./randist/binomial.c ./randist/binomial_tpe.c ./randist/cauchy.c ./randist/chisq.c ./randist/dirichlet.c ./randist/discrete.c ./randist/erlang.c ./randist/exponential.c ./randist/exppow.c ./randist/fdist.c ./randist/flat.c ./randist/gamma.c ./randist/gauss.o ./randist/gausstail.o ./randist/gausszig.o ./randist/geometric.o ./randist/gumbel.o ./randist/hyperg.o ./randist/landau.o ./randist/laplace.o ./randist/levy.o ./randist/logarithmic.o ./randist/logistic.o ./randist/lognormal.o ./randist/multinomial.o ./randist/mvgauss.o ./randist/nbinomial.o ./randist/pareto.o ./randist/pascal.o ./randist/poisson.o ./randist/rayleigh.o ./randist/shuffle.o ./randist/sphere.o ./randist/tdist.o ./randist/weibull.o ./sort/sort.o ./sort/sortind.o ./sort/sortvec.o ./sort/sortvecind.o ./sort/subset.o ./sort/subsetind.o  ./permutation/canonical.o ./permutation/file.o ./permutation/init.o ./permutation/inline.o ./permutation/permutation.o ./permutation/permute.o ./multimin/conjugate_fr.o ./multimin/conjugate_pr.o ./multimin/convergence.o ./multimin/diff.o ./multimin/fdfminimizer.o ./multimin/fminimizer.o ./multimin/linear_minimize.o ./multimin/simplex2.o ./multimin/simplex.o ./multimin/steepest_descent.o ./multimin/vector_bfgs2.o ./multimin/vector_bfgs.o  ./histogram/add2d.o ./histogram/add.o ./histogram/calloc_range2d.o ./histogram/calloc_range.o ./histogram/copy2d.o ./histogram/copy.o ./histogram/file2d.o ./histogram/file.o ./histogram/get2d.o ./histogram/get.o ./histogram/init2d.o ./histogram/init.o ./histogram/maxval2d.o ./histogram/maxval.o ./histogram/oper2d.o ./histogram/oper.o ./histogram/params2d.o ./histogram/params.o ./histogram/pdf2d.o ./histogram/pdf.o ./histogram/reset2d.o ./histogram/reset.o ./histogram/stat2d.o ./histogram/stat.o ./histogram/urand.o ./block/block.o ./block/file.o ./block/init.o  ./sys/coerce.o ./sys/expm1.o ./sys/fcmp.o ./sys/fdiv.o ./sys/hypot.o ./sys/infnan.o ./sys/invhyp.o ./sys/ldfrexp.o ./sys/log1p.o ./sys/minmax.o ./sys/pow_int.o ./sys/prec.o ./blas/blas.o


# Some files fail to compile for several reasons:
##./specfunc/cheb_eval.o ./specfunc/cheb_eval_mode.o ./specfunc/legendre_source.o ./specfunc/legendre_con.o ./specfunc/legendre_H3d.o ./specfunc/legendre_P.o ./specfunc/legendre_poly.o ./specfunc/legendre_Qn.o

all: 	gsl
	${AR} -q libda_gsl.a ${OBJECTS}
	${RANLIB} libda_gsl.a

gsl:
	cd ${GSL_VECTOR}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"; cd ../; cd ${GSL_MATRIX}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"; cd ../; cd ${GSL_ERR}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"; cd ../; cd ${GSL_RNG}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"; cd ../; cd ${GSL_RANDIST}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"; cd ../; cd ${GSL_CDF}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"; cd ../; cd ${GSL_SPECFUNC}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"; cd ../; cd ${GSL_SORT}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"; cd ../; cd ${GSL_PERMUTATION}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"; cd ../; cd ${GSL_MULTIMIN}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"; cd ../; cd ${GSL_HISTOGRAM}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"; cd ../; cd ${GSL_BLOCK}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"; cd ../; cd ${GSL_SYS}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"; cd ../; cd ${GSL_BLAS}; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS) $(CPICFLAGS)"

clean:
	cd ${GSL_VECTOR}; ${MAKE} clean; cd ../; cd ${GSL_MATRIX}; ${MAKE} clean; cd ../; cd ${GSL_ERR}; ${MAKE} clean; cd ../; cd ${GSL_RNG}; ${MAKE} clean; cd ../; cd ${GSL_RANDIST};  ${MAKE} clean; cd ../; cd ${GSL_CDF}; ${MAKE} clean; cd ../; cd ${GSL_SPECFUNC}; ${MAKE} clean; cd ../; cd ${GSL_SORT}; ${MAKE} clean; cd ../; cd ${GSL_PERMUTATION}; ${MAKE} clean; cd ../; cd ${GSL_MULTIMIN}; ${MAKE} clean; cd ../; cd ${GSL_HISTOGRAM}; ${MAKE} clean; cd ../; cd ${GSL_BLOCK}; ${MAKE} clean; cd ../; cd ${GSL_SYS}; ${MAKE} clean; cd ../; cd ${GSL_BLAS}; ${MAKE} clean
