// Copyright © 2017 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#include <vinecopulib/misc/tools_c.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

///////////////////////////////////////////////////////////////////////////////
//  Function that allocates space and creates a double matrix.
//  Input: Dimension of the matrix to be created//  Output: Pointer to the created matrix.
///////////////////////////////////////////////////////////////////////////////
double **create_matrix(int rows, int columns)
{
    double **a;
    int i=0;
    a = (double**) calloc(rows, sizeof(double*));
    for(i=0;i<rows;i++) a[i] = (double*) calloc(columns,sizeof(double));
    return a;
}

///////////////////////////////////////////////////////////////////////////////
//  Function that frees the space that a double matrix has been allocated.
//  Input: Dimension of the matrix and a pointer to the matrix.
//  Output: Void.
///////////////////////////////////////////////////////////////////////////////
void free_matrix(double **a, int rows)
{
    int i=0;
    for(i=0;i<rows;i++) free(a[i]);
    free(a);
}

////////////////////////////////////////////////////////////////////////
// Fast computation of Kendall's tau
//
// Given by Shing (Eric) Fu, Feng Zhu, Guang (Jack) Yang, and Harry Joe
// Based on work of the method by Knight (1966)
/////////////////////////////////////////////////////////////////////////

void ktau(double *X, double *Y, int *N, double *tau, double *S, double *D, int *T, int *U, int *V)
{
	// Defining variables
	int K, L, I, J, Iend, Jend;
	int i, j, m;
	double *Y2 = calloc(*N, sizeof(double));
	double *X2 = calloc(*N, sizeof(double));
	double *xptr,*yptr; // HJ addition for swapping
	bool Iflag, Jflag, Xflag;
	*S = 0.; *D = 0.; *T = 0; *U = 0; *V = 0;

	/* 1.1 Sort X and Y in X order */
	/* Break ties in X according to Y */
	K=1;
	do
	{
		L=0;
		do
		{
			I = L;
			J = (I+K)<(*N)?(I+K):(*N);
			Iend = J;
			Jend = (J+K)<(*N)?(J+K):(*N);
			do
			{
				Iflag = (I < Iend);
				Jflag = (J < Jend);
				if (Iflag & Jflag)
				{
				 	Xflag = ((X[I] > X[J]) | ((X[I] == X[J]) & (Y[I] > Y[J])));
				}
				else
				{
					Xflag = false;
				}
				if((Iflag & !Jflag) | (Iflag & Jflag & !Xflag))
				{
					X2[L] = X[I];
					Y2[L] = Y[I];
					I++;
					L++;
				};
				if((!Iflag && Jflag) | (Iflag && Jflag && Xflag))
				{
					X2[L] = X[J];
					Y2[L] = Y[J];
					J++;
					L++;
				};
			}
			while(Iflag | Jflag);
		}
		while(L < *N);

		// Swap lists
		xptr=X; X=X2; X2=xptr;
		yptr=Y; Y=Y2; Y2=yptr;
		#ifdef OLD
		for(i = 0; i < *N; i++)
		{
			Xtem = X[i]; Ytem = Y[i];
			X[i] = X2[i]; Y[i] = Y2[i];
			X2[i] = Xtem; Y2[i] = Ytem;
		};
		#endif
		K *= 2;
	}
	while (K < *N);

	/* 1.2 Count pairs of tied X, T */
	j = 1;
	m = 1;
	for(i = 1; i < *N; i++)
    if(X[i] == X[i-1])
    {
		j++;
		if(Y[i] == Y[i-1])
		m++;
    }
    else if(j > 1)
    {
      *T += j * (j - 1) / 2;
      if(m > 1)
		*V += m * (m - 1) / 2;
      j = 1;
      m = 1;
    };
	*T += j * (j - 1) / 2;
	*V += m * (m - 1) / 2;

	/* 2.1 Sort Y again and count exchanges, S */
	/* Keep original relative order if tied */

	K=1;
	do
	{
		L=0;
		do
		{
			I = L;
			J = (I+K)<(*N)?(I+K):(*N);
			Iend = J;
			Jend = (J+K)<(*N)?(J+K):(*N);
			do
			{
				Iflag = (I < Iend);
				Jflag = (J < Jend);
				if (Iflag & Jflag)
				{
				 	Xflag = (Y[I] > Y[J]);
				}
				else
				{
					Xflag = false;
				}
				if((Iflag & !Jflag) | (Iflag & Jflag & !Xflag))
				{
					X2[L] = X[I];
					Y2[L] = Y[I];
					I++;
					L++;
				};
				if((!Iflag && Jflag) | (Iflag && Jflag && Xflag))
				{
					X2[L] = X[J];
					Y2[L] = Y[J];
					*S += Iend - I;
					J++;
					L++;
				};
			}
			while((Iflag | Jflag));
		}
		while(L < *N);

		// Swap lists
		xptr=X; X=X2; X2=xptr;
		yptr=Y; Y=Y2; Y2=yptr;
		#ifdef OLD
		for(i = 0; i < *N; i++)
		{
			Xtem = X[i]; Ytem = Y[i];
			X[i] = X2[i]; Y[i] = Y2[i];
			X2[i] = Xtem; Y2[i] = Ytem;
		};
		#endif
		K *= 2;
	}
	while (K < *N);

	/* 2.2 Count pairs of tied Y, U */
	j=1;
	for(i = 1; i < *N; i++)
		if(Y[i] == Y[i-1])
			j++;
		else if(j > 1)
		{
			*U += j * (j - 1) / 2;
			j = 1;
		};
	*U += j * (j - 1) / 2;


	/* 3. Calc. Kendall's Score and Denominator */
	*D = 0.5 * (*N) * (*N - 1);
	*S = *D - (2. * (*S) + *T + *U - *V);
	//if(*T > 0 | *U > 0) // adjust for ties
    *D = sqrt((*D - *T) * (*D - *U));
	*tau = (*S) / (*D);


  free(Y2);
  free(X2);
}


//////////////////////////////////////////////////
// ktau_matrix
//
// Input:
// data			data vector
// d			data dimension 1
// N			data dimension 2
//
// Output:
// out			Kendall's tau Matrix (as vector)

void ktau_matrix_c(double *data, int *d, int *N, double *out)
{
	double **x, S=0.0, D=0.0, *X, *Y;
	int k=0, i, j, t, T=0, U=0, V=0;
	x = create_matrix(*d,*N);
	X = (double*) calloc(*N,sizeof(double));
	Y = (double*) calloc(*N,sizeof(double));

	for(i=0;i<*d;i++)
    {
		for (t=0;t<*N;t++ )
		{
			x[i][t] = data[k];
			k++;
		}
    }

	k=0;
	for(i=0;i<((*d)-1);i++)
	{
		for(j=(i+1);j<(*d);j++)
		{
			for(t=0;t<*N;t++)
			{
				X[t]=x[i][t];
				Y[t]=x[j][t];
			}
			ktau(X, Y, N, &out[k], &S, &D, &T, &U, &V);
			k++;
		}
	}

	free(X);free(Y);free_matrix(x, *d);
}

/** 1/(2pi) */
#define M_1_2PI .159154943091895335768883763373

/** Debye function of order n.
    int(t=0..x) t^n dt / [exp(t)-1]
 The underivative for n=0 is log[1-exp(-x)], which is infinite at x=0,
 so the corresponding Debye-function is not defined at n=0.

 Literature:
 Ng et al, Math. Comp. 24 (110) (1970) 405
 Guseinov et al, Intl. J. Thermophys. 28 (4) (2007) 1420
 Engeln et al, Colloid & Polymer Sci. 261 (9) (1983) 736
 Maximon , Proc. R. Soc. A 459 (2039) (2003) 2807

 @param[in] x the argument and upper limit of the integral. x>=0.
 @param[in] n the power in the numerator of the integral, 1<=n<=20 .
 @return the Debye function. Zero if x<=0, and -1 if n is outside the
     parameter range that is implemented.
 @author Richard J. Mathar
 @since 2007-10-31 implemented range n=8..10
*/
double debyen(const double x, const int n)
{
    if (x<=0. )
        return 0. ;
    if ( n <1 || n >20)
        return -1. ;

    /* for values up to 4.80 the list of zeta functions
    and the sum up to k < K are huge enough to gain
    numeric stability in the sum */
    if(x>= 3. )
    {
        double sum ;
        /* list of n! zeta(n+1) for n =0 up to the maximum n implemented.
        Limited by the cancellation of digits encountered at smaller x and larger n.
        Digits := 30 :
        for n from 1 to 30 do
                printf("%.29e, ", evalf(n!*Zeta(n+1))) ;
        od:
        */
        static double nzetan[] = { 0., 1.64493406684822643647241516665e+00, 2.40411380631918857079947632302e+00,
                                   6.49393940226682914909602217926e+00, 2.48862661234408782319527716750e+01,
                                   1.22081167438133896765742151575e+02, 7.26011479714984435324654235892e+02,
                                   5.06054987523763947046857360209e+03, 4.04009783987476348853278236554e+04,
                                   3.63240911422382626807143525567e+05, 3.63059331160662871299061884284e+06,
                                   3.99266229877310867023270732405e+07, 4.79060379889831452426876764501e+08,
                                   6.22740219341097176419285340896e+09, 8.71809578301720678451912203103e+10,
                                   1.30769435221891382089009990749e+12, 2.09229496794815109066316556880e+13,
                                   3.55688785859223715975612396717e+14, 6.40238592281892140073564945334e+15,
                                   1.21645216453639396669876696274e+17, 2.43290316850786132173725681824e+18,
                                   5.10909543543702856776502748606e+19, 1.12400086178089123060215294900e+21
        } ;

        /* constrained to the list of nzetan[] given above */
        if ( n >= sizeof(nzetan)/sizeof(double) )
            return -1. ;

        /* n!*zeta(n) is the integral for x=infinity , 27.1.3 */
        sum = nzetan[n] ;

        /* the number of terms needed in the k-sum for x=0,1,2,3...
        * Reflects the n=1 case, because higher n need less terms.
        */
        static int kLim[] = {0,0,0,13,10,8,7,6,5,5,4,4,4,3} ;

        const int kmax = ((int) x < sizeof(kLim)/sizeof(int)) ? kLim[(int)x] : 3 ;
        /* Abramowitz Stegun 27.1.2 */
        int k;
        for(k=1; k<=kmax ;k++)
        {
            /* do not use x(k+1)=xk+x to avoid loss of precision */
            const double xk = x*k ;
            double ksum= 1./xk ;
            double tmp = n*ksum/xk ;	/* n/(xk)^2 */
            int s;
            for (s=1 ; s <= n ; s++)
            {
                ksum += tmp ;
                tmp *= (n-s)/xk ;
            }
            sum -= exp(-xk)* ksum*pow(x,n+1.) ;
        }
        return sum ;
    }
    else
    {
        /* list of absolute values of Bernoulli numbers of index 2*k, multiplied  by
        (2*pi)^k/(2k)!, and 2 subtracted, k=0,1,2,3,4
        Digits := 60 :
        interface(prettyprint=0) :
        for k from 1 to 70 do
         printf("%.30e,\n",evalf( abs((2*Pi)^(2*k)*bernoulli(2*k)/(2*k)!)-2 )) ;
        od;
        */
        static double koeff[]= { 0., 1.289868133696452872944830333292e+00, 1.646464674222763830320073930823e-01,
                                 3.468612396889827942903585958184e-02, 8.154712395888678757370477017305e-03,
                                 1.989150255636170674291917800638e-03, 4.921731066160965972759960954793e-04,
                                 1.224962701174096585170902102707e-04, 3.056451881730374346514297527344e-05,
                                 7.634586529999679712923289243879e-06, 1.907924067745592226304077366899e-06,
                                 4.769010054554659800072963735060e-07, 1.192163781025189592248804158716e-07,
                                 2.980310965673008246931701326140e-08, 7.450668049576914109638408036805e-09,
                                 1.862654864839336365743529470042e-09, 4.656623667353010984002911951881e-10,
                                 1.164154417580540177848737197821e-10, 2.910384378208396847185926449064e-11,
                                 7.275959094757302380474472711747e-12, 1.818989568052777856506623677390e-12,
                                 4.547473691649305030453643155957e-13, 1.136868397525517121855436593505e-13,
                                 2.842170965606321353966861428348e-14, 7.105427382674227346596939068119e-15,
                                 1.776356842186163180619218277278e-15, 4.440892101596083967998640188409e-16,
                                 1.110223024969096248744747318102e-16, 2.775557561945046552567818981300e-17,
                                 6.938893904331845249488542992219e-18, 1.734723476023986745668411013469e-18,
                                 4.336808689994439570027820336642e-19, 1.084202172491329082183740080878e-19,
                                 2.710505431220232916297046799365e-20, 6.776263578041593636171406200902e-21,
                                 1.694065894509399669649398521836e-21, 4.235164736272389463688418879636e-22,
                                 1.058791184067974064762782460584e-22, 2.646977960169798160618902050189e-23,
                                 6.617444900424343177893912768629e-24, 1.654361225106068880734221123349e-24,
                                 4.135903062765153408791935838694e-25, 1.033975765691286264082026643327e-25,
                                 2.584939414228213340076225223666e-26, 6.462348535570530772269628236053e-27,
                                 1.615587133892632406631747637268e-27, 4.038967834731580698317525293132e-28,
                                 1.009741958682895139216954234507e-28, 2.524354896707237808750799932127e-29,
                                 6.310887241768094478219682436680e-30, 1.577721810442023614704107565240e-30,
                                 3.944304526105059031370476640000e-31, 9.860761315262647572437533499000e-32,
                                 2.465190328815661892443976898000e-32, 6.162975822039154730370601500000e-33,
                                 1.540743955509788682510501190000e-33, 3.851859888774471706184973900000e-34,
                                 9.629649721936179265360991000000e-35, 2.407412430484044816328953000000e-35,
                                 6.018531076210112040809600000000e-36, 1.504632769052528010200750000000e-36,
                                 3.761581922631320025497600000000e-37, 9.403954806578300063715000000000e-38,
                                 2.350988701644575015901000000000e-38, 5.877471754111437539470000000000e-39,
                                 1.469367938527859384580000000000e-39, 3.673419846319648458500000000000e-40,
                                 9.183549615799121117000000000000e-41, 2.295887403949780249000000000000e-41,
                                 5.739718509874450320000000000000e-42, 1.434929627468612270000000000000e-42
        } ;

        double sum=0. ;

        /* Abramowitz-Stegun 27.1.1 */
        const double x2pi=x*M_1_2PI ;
        int k;
        for(k=1;k< sizeof(koeff)/sizeof(double)-1 ;k++)
        {
            const double sumold=sum ;
            /* do not precompute x2pi^2 to avoid loss of precision */
            sum += (2.+koeff[k])*pow(x2pi,2.*k)/(2*k+n) ;
            k++ ;
            sum -= (2.+koeff[k])*pow(x2pi,2.*k)/(2*k+n) ;
            if(sum == sumold)
                break ;
        }
        sum += 1./n-x/(2*(1+n)) ;
        return sum*pow(x,(double)n) ;
    }
}

#undef M_1_2PI

// borrowed from f2c.h
typedef int (*U_fp)();
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
double d_sign(const double *a, const double *b)
{
	double x = fabs(*a);
	return *b >= 0 ? x : -x;
}

static int c__1 = 1;
static int c_true = 1;
static int c__0 = 0;
static double c_b91 = 1.;
static double c_b99 = 2.;
static int c__100 = 100;

/* A subroutine for computing non-central multivariate t probabilities. */
/* This subroutine uses an algorithm (QRSVN) described in the paper */
/* "Comparison of Methods for the Computation of Multivariate */
/* t-Probabilities", by Alan Genz and Frank Bretz */
/* J. Comp. Graph. Stat. 11 (2002), pp. 950-971. */

/* Parameters */

/* N      INTEGER, the number of variables. */
/* NU     INTEGER, the number of degrees of freedom. */
/* If NU < 1, then an MVN probability is computed. */
/* LOWER  DOUBLE PRECISION, array of lower integration limits. */
/* UPPER  DOUBLE PRECISION, array of upper integration limits. */
/* INFIN  INTEGER, array of integration limits ints: */
/* if INFIN(I) < 0, Ith limits are (-infinity, infinity); */
/* if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; */
/* if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); */
/* if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. */
/* CORREL DOUBLE PRECISION, array of correlation coefficients; */
/* the correlation coefficient in row I column J of the */
/* correlation matrixshould be stored in */
/* CORREL( J + ((I-2)*(I-1))/2 ), for J < I. */
/* The correlation matrix must be positive semi-definite. */
/* DELTA  DOUBLE PRECISION, array of non-centrality parameters. */
/* MAXPTS INTEGER, maximum number of function values allowed. This */
/* parameter can be used to limit the time. A sensible */
/* strategy is to start with MAXPTS = 1000*N, and then */
/* increase MAXPTS if ERROR is too large. */
/* ABSEPS DOUBLE PRECISION absolute error tolerance. */
/* RELEPS DOUBLE PRECISION relative error tolerance. */
/* ERROR  DOUBLE PRECISION estimated absolute error, */
/* with 99% confidence level. */
/* VALUE  DOUBLE PRECISION estimated value for the integral */
/* INFORM INTEGER, termination status parameter: */
/* if INFORM = 0, normal completion with ERROR < EPS; */
/* if INFORM = 1, completion with ERROR > EPS and MAXPTS */
/*           function vaules used; increase MAXPTS to */
/*           decrease ERROR; */
/* if INFORM = 2, N > 1000 or N < 1. */
/* if INFORM = 3, correlation matrix not positive semi-definite. */
void mvtdst_(int *n, int *nu, double *lower,
			double *upper, int *infin, double *correl, double *delta,
			int *maxpts, double *abseps, double *releps,
			double *error, double *value, int *inform__)
{
	static double e[1], v[1];
	static int nd;
	extern int mvkbrv_(int *, int *, int *, int *,
					   int (*)(), double *, double *, double *,
					   double *, int *);
	extern int mvsubr_();
	extern int mvints_(int *, int *, double *, double *,
					   double *, double *, int *, int *,
					   double *, double *, int *);

	/* Parameter adjustments */
	--delta;
	--correl;
	--infin;
	--upper;
	--lower;

	/* Function Body */
	int ivls = 0;
	if (*n > 1000 || *n < 1) {
		*value = 0.;
		*error = 1.;
		*inform__ = 2;
	} else {
		mvints_(n, nu, &correl[1], &lower[1], &upper[1], &delta[1], &infin[1],
				&nd, value, error, inform__);
		if (*inform__ == 0 && nd > 0) {

			/* Call the lattice rule integration subroutine */
			mvkbrv_(&nd, &ivls, maxpts, &c__1, (U_fp)mvsubr_, abseps,
					releps, e, v, inform__);
			*error = e[0];
			*value = v[0];
		}
	}
}


int mvsubr_0_(int n__, int *n, double *w, int *
nf, double *f, int *nuin, double *correl, double *
lower, double *upper, double *delta, int *infin, int *
nd, double *vl, double *er, int *inform__)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	static double a[1000], b[1000], r__, y[1000], di, ei, dl[1000];
	static int nu, ny;
	static double cov[500500], snu;
	static int infi[1000];
	extern double mvchnv_(int *, double *);
	extern int mvspcl_(int *, int *, double *,
					   double *, double *, double *, int *, double *,
					   double *, double *, int *);
	extern int mvvlsb_(int *,double *, double *, double *, int *,
					   double *, double *, double *, double *, double *,
					   double *, int *, double *);
	extern int mvsort_(int *, double *, double *, double *, double *,
					   int *, double *, int *, int *, double *,
					   double *, double *, double *, int *, int *);

	/* Integrand subroutine */
	/* Parameter adjustments */
	if (w) {
		--w;
	}
	if (f) {
		--f;
	}
	if (correl) {
		--correl;
	}
	if (lower) {
		--lower;
	}
	if (upper) {
		--upper;
	}
	if (delta) {
		--delta;
	}
	if (infin) {
		--infin;
	}

	/* Function Body */
	switch(n__) {
		case 1: goto L_mvints;
	}

	if (nu <= 0) {
		r__ = 1.;
		i__1 = *n + 1;
		mvvlsb_(&i__1, &w[1], &r__, dl, infi, a, b, cov, y, &di, &ei, &ny, &f[
				1]);
	} else {
		r__ = mvchnv_(&nu, &w[*n]) / snu;
		mvvlsb_(n, &w[1], &r__, dl, infi, a, b, cov, y, &di, &ei, &ny, &f[1]);
	}
	return 0;

	/* Entry point for intialization. */
	L_mvints:

	/* Initialization and computation of covariance Cholesky factor. */
	mvsort_(n, &lower[1], &upper[1], &delta[1], &correl[1], &infin[1], y, &
			c_true, nd, a, b, dl, cov, infi, inform__);
	nu = *nuin;
	mvspcl_(nd, &nu, a, b, dl, cov, infi, &snu, vl, er, inform__);
	return 0;
}

int mvsubr_(int *n, double *w, int *nf,
			double *f)
{
	return mvsubr_0_(0, n, w, nf, f, (int *)0, (double *)0, (double *)0,
					 (double *)0, (double *)0, (int *)0, (int *)0,
					 (double *)0, (double *)0, (int *)0);
}

int mvints_(int *n, int *nuin, double *correl,
			double *lower, double *upper, double *delta, int *
infin, int *nd, double *vl, double *er, int *inform__)
{
	return mvsubr_0_(1, n, (double *)0, (int *)0, (double *)0, nuin,
					 correl, lower, upper, delta, infin, nd, vl, er, inform__);
}


int mvspcl_(int *nd, int *nu, double *a,
			double *b, double *dl, double *cov, int *infi,
			double *snu, double *vl, double *er, int *inform__)
{
	/* System generated locals */
	double d__1;

	/* Builtin functions */
	double sqrt(double);

	/* Local variables */
	static double r__;
	extern double mvbvt_(int *, double *, double *, int *,
						 double *), mvstdt_(int *, double *);


	/* Special cases subroutine */
	/* Parameter adjustments */
	--infi;
	--cov;
	--dl;
	--b;
	--a;

	/* Function Body */
	if (*inform__ > 0) {
		*vl = 0.;
		*er = 1.;
	} else {

		/* Special cases */
		if (*nd == 0) {
			*er = 0.;
			/* Code added to fix ND = 0 bug, 24/03/2009 -> */
			*vl = 1.;
			/* <- Code added to fix ND = 0 bug, 24/03/2009 */
		} else if (*nd == 1 && (*nu < 1 || fabs(dl[1]) == 0.)) {

			/* 1-d case for normal or central t */
			*vl = 1.;
			if (infi[1] != 1) {
				d__1 = b[1] - dl[1];
				*vl = mvstdt_(nu, &d__1);
			}
			if (infi[1] != 0) {
				d__1 = a[1] - dl[1];
				*vl -= mvstdt_(nu, &d__1);
			}
			if (*vl < 0.) {
				*vl = 0.;
			}
			*er = 2e-16;
			*nd = 0;
		} else if (*nd == 2 && (*nu < 1 || fabs(dl[1]) + fabs(dl[2]) == 0.)) {

			/* 2-d case for normal or central t */
			if (infi[1] != 0) {
				a[1] -= dl[1];
			}
			if (infi[1] != 1) {
				b[1] -= dl[1];
			}
			if (infi[2] != 0) {
				a[2] -= dl[2];
			}
			if (infi[2] != 1) {
				b[2] -= dl[2];
			}
			if (fabs(cov[3]) > 0.) {

				/* 2-d nonsingular case */

				/* Computing 2nd power */
				d__1 = cov[2];
				r__ = sqrt(d__1 * d__1 + 1);
				if (infi[2] != 0) {
					a[2] /= r__;
				}
				if (infi[2] != 1) {
					b[2] /= r__;
				}
				cov[2] /= r__;
				*vl = mvbvt_(nu, &a[1], &b[1], &infi[1], &cov[2]);
				*er = 1e-15;
			} else {

				/* 2-d singular case */
				if (infi[1] != 0) {
					if (infi[2] != 0) {
						a[1] = max(a[1],a[2]);
					}
				} else {
					if (infi[2] != 0) {
						a[1] = a[2];
					}
				}
				if (infi[1] != 1) {
					if (infi[2] != 1) {
						b[1] = min(b[1],b[2]);
					}
				} else {
					if (infi[2] != 1) {
						b[1] = b[2];
					}
				}
				if (infi[1] != infi[2]) {
					infi[1] = 2;
				}
				*vl = 1.;
				/* A(1), B(1) Bug Fixed, 28/05/2013 */
				if (infi[1] != 1) {
					*vl = mvstdt_(nu, &b[1]);
				}
				if (infi[1] != 0) {
					*vl -= mvstdt_(nu, &a[1]);
				}
				if (*vl < 0.) {
					*vl = 0.;
				}
				*er = 2e-16;
			}
			*nd = 0;
		} else {
			if (*nu > 0) {
				*snu = sqrt((double) (*nu));
			} else {
				--(*nd);
			}
		}
	}
	return 0;
}

int mvvlsb_(int *n, double *w, double *r__,
			double *dl, int *infi, double *a, double *b,
			double *cov, double *y, double *di, double *ei,
			int *nd, double *value)
{
	/* System generated locals */
	int i__1, i__2;
	double d__1, d__2;

	/* Local variables */
	static int i__, j;
	static double ai, bi;
	static int ij;
	static double sum;
	static int infa, infb;
	extern int mvlims_(double *, double *, int *,
					   double *, double *);
	extern double mvphnv_(double *);


	/* Integrand subroutine */
	/* Parameter adjustments */
	--y;
	--cov;
	--b;
	--a;
	--infi;
	--dl;
	--w;

	/* Function Body */
	*value = 1.;
	infa = 0;
	infb = 0;
	*nd = 0;
	ij = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = dl[i__];
		i__2 = i__ - 1;
		for (j = 1; j <= i__2; ++j) {
			++ij;
			if (j <= *nd) {
				sum += cov[ij] * y[j];
			}
		}
		if (infi[i__] != 0) {
			if (infa == 1) {
				/* Computing MAX */
				d__1 = ai, d__2 = *r__ * a[i__] - sum;
				ai = max(d__1,d__2);
			} else {
				ai = *r__ * a[i__] - sum;
				infa = 1;
			}
		}
		if (infi[i__] != 1) {
			if (infb == 1) {
				/* Computing MIN */
				d__1 = bi, d__2 = *r__ * b[i__] - sum;
				bi = min(d__1,d__2);
			} else {
				bi = *r__ * b[i__] - sum;
				infb = 1;
			}
		}
		++ij;
		if (i__ == *n || cov[ij + *nd + 2] > 0.) {
			i__2 = infa + infa + infb - 1;
			mvlims_(&ai, &bi, &i__2, di, ei);
			if (*di >= *ei) {
				*value = 0.;
				return 0;
			} else {
				*value *= *ei - *di;
				++(*nd);
				if (i__ < *n) {
					d__1 = *di + w[*nd] * (*ei - *di);
					y[*nd] = mvphnv_(&d__1);
				}
				infa = 0;
				infb = 0;
			}
		}
	}
	return 0;
}

/* Subroutine to sort integration limits and determine Cholesky factor. */
int mvsort_(int *n, double *lower, double *upper,
			double *delta, double *correl, int *infin, double *y,
			int *pivot, int *nd, double *a, double *b,
			double *dl, double *cov, int *infi, int *inform__)
{
	/* System generated locals */
	int i__1, i__2, i__3, i__4;
	double d__1;

	/* Builtin functions */
	double sqrt(double);

	/* Local variables */
	static double d__, e;
	static int i__, j, k, l, m;
	static double aj, bj;
	static int ii, ij, il, jl;
	static double sum, amin, bmin;
	static int jmin;
	static double epsi, demin, sumsq, cvdiag;
	extern int mvlims_(double *, double *, int *,
					   double *, double *);
	extern double mvtdns_(int *, double *);
	extern int mvswap_(int *, int *, double *, double *, double *,
					   int *, int *, double *),
			mvsswp_(double *, double *);

	/* Parameter adjustments */
	--infi;
	--cov;
	--dl;
	--b;
	--a;
	--y;
	--infin;
	--correl;
	--delta;
	--upper;
	--lower;

	/* Function Body */
	*inform__ = 0;
	ij = 0;
	ii = 0;
	*nd = *n;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		a[i__] = 0.;
		b[i__] = 0.;
		dl[i__] = 0.;
		infi[i__] = infin[i__];
		if (infi[i__] < 0) {
			--(*nd);
		} else {
			if (infi[i__] != 0) {
				a[i__] = lower[i__];
			}
			if (infi[i__] != 1) {
				b[i__] = upper[i__];
			}
			dl[i__] = delta[i__];
		}
		i__2 = i__ - 1;
		for (j = 1; j <= i__2; ++j) {
			++ij;
			++ii;
			cov[ij] = correl[ii];
		}
		++ij;
		cov[ij] = 1.;
	}

	/* First move any doubly infinite limits to innermost positions. */
	if (*nd > 0) {
		i__1 = *nd + 1;
		for (i__ = *n; i__ >= i__1; --i__) {
			if (infi[i__] >= 0) {
				i__2 = i__ - 1;
				for (j = 1; j <= i__2; ++j) {
					if (infi[j] < 0) {
						mvswap_(&j, &i__, &a[1], &b[1], &dl[1], &infi[1], n, &
								cov[1]);
						goto L10;
					}
				}
			}
			L10:
			;
		}

		/* Sort remaining limits and determine Cholesky factor. */
		ii = 0;
		jl = *nd;
		i__1 = *nd;
		for (i__ = 1; i__ <= i__1; ++i__) {

			/* Determine the integration limits for variable with minimum */
			/* expected probability and interchange that variable with Ith. */
			demin = 1.;
			jmin = i__;
			cvdiag = 0.;
			ij = ii;
			epsi = i__ * 1e-10;
			if (! (*pivot)) {
				jl = i__;
			}
			i__2 = jl;
			for (j = i__; j <= i__2; ++j) {
				if (cov[ij + j] > epsi) {
					sumsq = sqrt(cov[ij + j]);
					sum = dl[j];
					i__3 = i__ - 1;
					for (k = 1; k <= i__3; ++k) {
						sum += cov[ij + k] * y[k];
					}
					aj = (a[j] - sum) / sumsq;
					bj = (b[j] - sum) / sumsq;
					mvlims_(&aj, &bj, &infi[j], &d__, &e);
					if (demin >= e - d__) {
						jmin = j;
						amin = aj;
						bmin = bj;
						demin = e - d__;
						cvdiag = sumsq;
					}
				}
				ij += j;
			}
			if (jmin > i__) {
				mvswap_(&i__, &jmin, &a[1], &b[1], &dl[1], &infi[1], n, &cov[
						1]);
			}
			if (cov[ii + i__] < -epsi) {
				*inform__ = 3;
			}
			cov[ii + i__] = cvdiag;

			/* Compute Ith column of Cholesky factor. */
			/* Compute expected value for Ith integration variable and */
			/* scale Ith covariance matrix row and limits. */
			if (cvdiag > 0.) {
				il = ii + i__;
				i__2 = *nd;
				for (l = i__ + 1; l <= i__2; ++l) {
					cov[il + i__] /= cvdiag;
					ij = ii + i__;
					i__3 = l;
					for (j = i__ + 1; j <= i__3; ++j) {
						cov[il + j] -= cov[il + i__] * cov[ij + i__];
						ij += j;
					}
					il += l;
				}

				/* Expected Y = -( density(b) - density(a) )/( b - a ) */
				if (demin > epsi) {
					y[i__] = 0.;
					if (infi[i__] != 0) {
						y[i__] = mvtdns_(&c__0, &amin);
					}
					if (infi[i__] != 1) {
						y[i__] -= mvtdns_(&c__0, &bmin);
					}
					y[i__] /= demin;
				} else {
					if (infi[i__] == 0) {
						y[i__] = bmin;
					}
					if (infi[i__] == 1) {
						y[i__] = amin;
					}
					if (infi[i__] == 2) {
						y[i__] = (amin + bmin) / 2;
					}
				}
				i__2 = i__;
				for (j = 1; j <= i__2; ++j) {
					++ii;
					cov[ii] /= cvdiag;
				}
				a[i__] /= cvdiag;
				b[i__] /= cvdiag;
				dl[i__] /= cvdiag;
			} else {
				il = ii + i__;
				i__2 = *nd;
				for (l = i__ + 1; l <= i__2; ++l) {
					cov[il + i__] = 0.;
					il += l;
				}

				/* If the covariance matrix diagonal entry is zero, */
				/* permute limits and rows, if necessary. */
				for (j = i__ - 1; j >= 1; --j) {
					if ((d__1 = cov[ii + j], fabs(d__1)) > epsi) {
						a[i__] /= cov[ii + j];
						b[i__] /= cov[ii + j];
						dl[i__] /= cov[ii + j];
						if (cov[ii + j] < 0.) {
							mvsswp_(&a[i__], &b[i__]);
							if (infi[i__] != 2) {
								infi[i__] = 1 - infi[i__];
							}
						}
						i__2 = j;
						for (l = 1; l <= i__2; ++l) {
							cov[ii + l] /= cov[ii + j];
						}
						i__2 = i__ - 1;
						for (l = j + 1; l <= i__2; ++l) {
							if (cov[(l - 1) * l / 2 + j + 1] > 0.) {
								ij = ii;
								i__3 = l;
								for (k = i__ - 1; k >= i__3; --k) {
									i__4 = k;
									for (m = 1; m <= i__4; ++m) {
										mvsswp_(&cov[ij - k + m], &cov[ij + m]
										);
									}
									mvsswp_(&a[k], &a[k + 1]);
									mvsswp_(&b[k], &b[k + 1]);
									mvsswp_(&dl[k], &dl[k + 1]);
									m = infi[k];
									infi[k] = infi[k + 1];
									infi[k + 1] = m;
									ij -= k;
								}
								goto L20;
							}
						}
						goto L20;
					}
					cov[ii + j] = 0.;
				}
				L20:
				ii += i__;
				y[i__] = 0.;
			}
		}
	}
	return 0;
}


double mvtdns_(int *nu, double *x)
{
	/* System generated locals */
	int i__1;
	double ret_val, d__1;

	/* Builtin functions */
	double sqrt(double), exp(double);

	/* Local variables */
	static int i__;
	static double prod;

	ret_val = 0.;
	if (*nu > 0) {
		prod = 1 / sqrt((double) (*nu));
		for (i__ = *nu - 2; i__ >= 1; i__ += -2) {
			prod = prod * (i__ + 1) / i__;
		}
		if (*nu % 2 == 0) {
			prod /= 2;
		} else {
			prod /= 3.141592653589793;
		}
		d__1 = sqrt(*x * *x / *nu + 1);
		i__1 = *nu + 1;
		ret_val = prod / pow(d__1, i__1);
	} else {
		if (fabs(*x) < 10.) {
			ret_val = exp(-(*x) * *x / 2) / 2.506628274631001;
		}
	}
	return ret_val;
}


int mvlims_(double *a, double *b, int *infin,
			double *lower, double *upper)
{
	extern double mvphi_(double *);

	*lower = 0.;
	*upper = 1.;
	if (*infin >= 0) {
		if (*infin != 0) {
			*lower = mvphi_(a);
		}
		if (*infin != 1) {
			*upper = mvphi_(b);
		}
	}
	*upper = max(*upper,*lower);
	return 0;
}


int mvsswp_(double *x, double *y)
{
	static double t;

	t = *x;
	*x = *y;
	*y = t;
	return 0;
}

/* Swaps rows and columns P and Q in situ, with P <= Q. */
int mvswap_(int *p, int *q, double *a,
			double *b, double *d__, int *infin, int *n,
			double *c__)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	static int i__, j, ii, jj;
	extern int mvsswp_(double *, double *);

	/* Parameter adjustments */
	--c__;
	--infin;
	--d__;
	--b;
	--a;

	/* Function Body */
	mvsswp_(&a[*p], &a[*q]);
	mvsswp_(&b[*p], &b[*q]);
	mvsswp_(&d__[*p], &d__[*q]);
	j = infin[*p];
	infin[*p] = infin[*q];
	infin[*q] = j;
	jj = *p * (*p - 1) / 2;
	ii = *q * (*q - 1) / 2;
	mvsswp_(&c__[jj + *p], &c__[ii + *q]);
	i__1 = *p - 1;
	for (j = 1; j <= i__1; ++j) {
		mvsswp_(&c__[jj + j], &c__[ii + j]);
	}
	jj += *p;
	i__1 = *q - 1;
	for (i__ = *p + 1; i__ <= i__1; ++i__) {
		mvsswp_(&c__[jj + *p], &c__[ii + i__]);
		jj += i__;
	}
	ii += *q;
	i__1 = *n;
	for (i__ = *q + 1; i__ <= i__1; ++i__) {
		mvsswp_(&c__[ii + *p], &c__[ii + *q]);
		ii += i__;
	}
	return 0;
}

/* Normal distribution probabilities accurate to 1d-15. */
/* Reference: J.L. Schonfelder, Math Comp 32(1978), pp 1232-1240. */
double mvphi_(double *z__)
{
	/* Initialized data */

	static double a[44] = { .610143081923200417926465815756,
							-.434841272712577471828182820888,.176351193643605501125840298123,
							-.060710795609249414860051215825,.017712068995694114486147141191,
							-.004321119385567293818599864968,8.54216676887098678819832055e-4,
							-1.2715509060916274262889394e-4,1.1248167243671189468847072e-5,
							3.13063885421820972630152e-7,-2.70988068537762022009086e-7,
							3.0737622701407688440959e-8,2.515620384817622937314e-9,
							-1.02892992132031912759e-9,2.9944052119949939363e-11,
							2.605178968726693629e-11,-2.634839924171969386e-12,
							-6.43404509890636443e-13,1.12457401801663447e-13,
							1.7281533389986098e-14,-4.264101694942375e-15,
							-5.45371977880191e-16,1.58697607761671e-16,2.0899837844334e-17,
							-5.900526869409e-18,-9.41893387554e-19,2.1497735647e-19,
							4.6660985008e-20,-7.243011862e-21,-2.387966824e-21,1.91177535e-22,
							1.20482568e-22,-6.72377e-25,-5.747997e-24,-4.28493e-25,
							2.44856e-25,4.3793e-26,-8.151e-27,-3.089e-27,9.3e-29,1.74e-28,
							1.6e-29,-8e-30,-2e-30 };

	/* System generated locals */
	double ret_val;

	/* Builtin functions */
	double exp(double);

	/* Local variables */
	static double b;
	static int i__;
	static double p, t, bm, bp, xa;

	xa = fabs(*z__) / 1.414213562373095048801688724209;
	if (xa > 100.) {
		p = 0.;
	} else {
		t = (xa * 8 - 30) / (xa * 4 + 15);
		bm = 0.;
		b = 0.;
		for (i__ = 24; i__ >= 0; --i__) {
			bp = b;
			b = bm;
			bm = t * b - bp + a[i__];
		}
		p = exp(-xa * xa) * (bm - bp) / 4;
	}
	if (*z__ > 0.) {
		p = 1 - p;
	}
	ret_val = p;
	return ret_val;
}

/* 	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3 */

/* 	Produces the normal deviate Z corresponding to a given lower */
/* 	tail area of P. */
/* 	The hash sums below are the sums of the mantissas of the */
/* 	coefficients.   They are included for use in checking */
/* 	transcription. */
/* Coefficients for P close to 0.5 */
/* HASH SUM AB    55.88319 28806 14901 4439 */
/* Coefficients for P not close to 0, 0.5 or 1. */
/* HASH SUM CD    49.33206 50330 16102 89036 */
/* 	   Coefficients for P near 0 or 1. */
/* HASH SUM EF    47.52583 31754 92896 71629 */
double mvphnv_(double *p)
{
	/* System generated locals */
	double ret_val, d__1, d__2;

	/* Builtin functions */
	double log(double), sqrt(double);

	/* Local variables */
	static double q, r__;

	q = (*p * 2 - 1) / 2;
	if (fabs(q) <= .425) {
		r__ = .180625 - q * q;
		ret_val = q * (((((((r__ * 2509.0809287301226727 +
							 33430.575583588128105) * r__ + 67265.770927008700853) * r__ +
						   45921.953931549871457) * r__ + 13731.693765509461125) * r__ +
						 1971.5909503065514427) * r__ + 133.14166789178437745) * r__ +
					   3.387132872796366608) / (((((((r__ * 5226.495278852854561 +
													  28729.085735721942674) * r__ + 39307.89580009271061) * r__ +
													21213.794301586595867) * r__ + 5394.1960214247511077) * r__ +
												  687.1870074920579083) * r__ + 42.313330701600911252) * r__ +
												1);
	} else {
		/* Computing MIN */
		d__1 = *p, d__2 = 1 - *p;
		r__ = min(d__1,d__2);
		if (r__ > 0.) {
			r__ = sqrt(-log(r__));
			if (r__ <= 5.) {
				r__ += -1.6;
				ret_val = (((((((r__ * 7.7454501427834140764e-4 +
								 .0227238449892691845833) * r__ +
								.24178072517745061177) * r__ + 1.27045825245236838258)
							  * r__ + 3.64784832476320460504) * r__ +
							 5.7694972214606914055) * r__ + 4.6303378461565452959)
						   * r__ + 1.42343711074968357734) / (((((((r__ *
																	1.05075007164441684324e-9 + 5.475938084995344946e-4) *
																   r__ + .0151986665636164571966) * r__ +
																  .14810397642748007459) * r__ + .68976733498510000455)
																* r__ + 1.6763848301838038494) * r__ +
															   2.05319162663775882187) * r__ + 1);
			} else {
				r__ += -5.;
				ret_val = (((((((r__ * 2.01033439929228813265e-7 +
								 2.71155556874348757815e-5) * r__ +
								.0012426609473880784386) * r__ +
							   .026532189526576123093) * r__ + .29656057182850489123)
							 * r__ + 1.7848265399172913358) * r__ +
							5.4637849111641143699) * r__ + 6.6579046435011037772)
						  / (((((((r__ * 2.04426310338993978564e-15 +
								   1.4215117583164458887e-7) * r__ +
								  1.8463183175100546818e-5) * r__ +
								 7.868691311456132591e-4) * r__ +
								.0148753612908506148525) * r__ +
							   .13692988092273580531) * r__ + .59983220655588793769)
							 * r__ + 1);
			}
		} else {
			ret_val = 9.;
		}
		if (q < 0.) {
			ret_val = -ret_val;
		}
	}
	return ret_val;
}

/* A function for computing bivariate normal probabilities. */

/* Parameters */

/* LOWER  REAL, array of lower integration limits. */
/* UPPER  REAL, array of upper integration limits. */
/* INFIN  INTEGER, array of integration limits ints: */
/* if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; */
/* if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); */
/* if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. */
/* CORREL REAL, correlation coefficient. */
double mvbvn_(double *lower, double *upper, int *infin,
			  double *correl)
{
	/* System generated locals */
	double ret_val, d__1, d__2, d__3, d__4;

	/* Local variables */
	extern double mvbvu_(double *, double *, double *);

	/* Parameter adjustments */
	--infin;
	--upper;
	--lower;

	/* Function Body */
	if (infin[1] == 2 && infin[2] == 2) {
		ret_val = mvbvu_(&lower[1], &lower[2], correl) - mvbvu_(&upper[1], &
				lower[2], correl) - mvbvu_(&lower[1], &upper[2], correl) +
				  mvbvu_(&upper[1], &upper[2], correl);
	} else if (infin[1] == 2 && infin[2] == 1) {
		ret_val = mvbvu_(&lower[1], &lower[2], correl) - mvbvu_(&upper[1], &
				lower[2], correl);
	} else if (infin[1] == 1 && infin[2] == 2) {
		ret_val = mvbvu_(&lower[1], &lower[2], correl) - mvbvu_(&lower[1], &
				upper[2], correl);
	} else if (infin[1] == 2 && infin[2] == 0) {
		d__1 = -upper[1];
		d__2 = -upper[2];
		d__3 = -lower[1];
		d__4 = -upper[2];
		ret_val = mvbvu_(&d__1, &d__2, correl) - mvbvu_(&d__3, &d__4, correl);
	} else if (infin[1] == 0 && infin[2] == 2) {
		d__1 = -upper[1];
		d__2 = -upper[2];
		d__3 = -upper[1];
		d__4 = -lower[2];
		ret_val = mvbvu_(&d__1, &d__2, correl) - mvbvu_(&d__3, &d__4, correl);
	} else if (infin[1] == 1 && infin[2] == 0) {
		d__1 = -upper[2];
		d__2 = -(*correl);
		ret_val = mvbvu_(&lower[1], &d__1, &d__2);
	} else if (infin[1] == 0 && infin[2] == 1) {
		d__1 = -upper[1];
		d__2 = -(*correl);
		ret_val = mvbvu_(&d__1, &lower[2], &d__2);
	} else if (infin[1] == 1 && infin[2] == 1) {
		ret_val = mvbvu_(&lower[1], &lower[2], correl);
	} else if (infin[1] == 0 && infin[2] == 0) {
		d__1 = -upper[1];
		d__2 = -upper[2];
		ret_val = mvbvu_(&d__1, &d__2, correl);
	} else {
		ret_val = 1.;
	}
	return ret_val;
}

/* A function for computing bivariate normal probabilities; */
/* developed using */
/* Drezner, Z. and Wesolowsky, G. O. (1989), */
/* On the Computation of the Bivariate Normal Integral, */
/* J. Stat. Comput. Simul.. 35 pp. 101-107. */
/* with extensive modications for double precisions by */
/* Alan Genz and Yihong Ge */
/* Department of Mathematics */
/* Washington State University */
/* Pullman, WA 99164-3113 */
/* Email : alangenz@wsu.edu */

/* BVN - calculate the probability that X is larger than SH and Y is */
/* larger than SK. */

/* Parameters */

/* SH  REAL, integration limit */
/* SK  REAL, integration limit */
/* R   REAL, correlation coefficient */
/* LG  INTEGER, number of Gauss Rule Points and Weights */

/* Gauss Legendre Points and Weights, N =  6 */
/* Gauss Legendre Points and Weights, N = 12 */
/* Gauss Legendre Points and Weights, N = 20 */
double mvbvu_(double *sh, double *sk, double *r__)
{
	/* Initialized data */

	static struct {
		double e_1[3];
		double fill_2[7];
		double e_3[6];
		double fill_4[4];
		double e_5[10];
	} equiv_112 = { .1713244923791705, .3607615730481384,
					.4679139345726904, {0}, .04717533638651177, .1069393259953183,
					.1600783285433464, .2031674267230659, .2334925365383547,
					.2491470458134029, {0}, .01761400713915212,
					.04060142980038694, .06267204833410906, .08327674157670475,
					.1019301198172404, .1181945319615184, .1316886384491766,
					.1420961093183821, .1491729864726037, .1527533871307259 };

#define w ((double *)&equiv_112)

	static struct {
		double e_1[3];
		double fill_2[7];
		double e_3[6];
		double fill_4[4];
		double e_5[10];
	} equiv_113 = { -.9324695142031522, -.6612093864662647,
					-.238619186083197, {0}, -.9815606342467191, -.904117256370475,
					-.769902674194305, -.5873179542866171, -.3678314989981802,
					-.1252334085114692, {0}, -.9931285991850949,
					-.9639719272779138, -.9122344282513259, -.8391169718222188,
					-.7463319064601508, -.636053680726515, -.5108670019508271,
					-.3737060887154196, -.2277858511416451, -.07652652113349733 };

#define x ((double *)&equiv_113)


	/* System generated locals */
	int i__1;
	double ret_val, d__1, d__2;

	/* Builtin functions */
	double asin(double), sin(double), exp(double), sqrt(
			double);

	/* Local variables */
	static double a, b, c__, d__, h__;
	static int i__;
	static double k;
	static int lg;
	static double as;
	static int ng;
	static double bs, hk, hs, sn, rs, xs, bvn, asr;
	extern double mvphi_(double *);


	if (fabs(*r__) < .3f) {
		ng = 1;
		lg = 3;
	} else if (fabs(*r__) < .75f) {
		ng = 2;
		lg = 6;
	} else {
		ng = 3;
		lg = 10;
	}
	h__ = *sh;
	k = *sk;
	hk = h__ * k;
	bvn = 0.;
	if (fabs(*r__) < .925f) {
		hs = (h__ * h__ + k * k) / 2;
		asr = asin(*r__);
		i__1 = lg;
		for (i__ = 1; i__ <= i__1; ++i__) {
			sn = sin(asr * (x[i__ + ng * 10 - 11] + 1) / 2);
			bvn += w[i__ + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn))
					;
			sn = sin(asr * (-x[i__ + ng * 10 - 11] + 1) / 2);
			bvn += w[i__ + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn))
					;
		}
		d__1 = -h__;
		d__2 = -k;
		bvn = bvn * asr / 12.566370614359172 + mvphi_(&d__1) * mvphi_(&d__2);
	} else {
		if (*r__ < 0.) {
			k = -k;
			hk = -hk;
		}
		if (fabs(*r__) < 1.) {
			as = (1 - *r__) * (*r__ + 1);
			a = sqrt(as);
			/* Computing 2nd power */
			d__1 = h__ - k;
			bs = d__1 * d__1;
			c__ = (4 - hk) / 8;
			d__ = (12 - hk) / 16;
			bvn = a * exp(-(bs / as + hk) / 2) * (1 - c__ * (bs - as) * (1 -
																		 d__ * bs / 5) / 3 + c__ * d__ * as * as / 5);
			if (hk > -160.) {
				b = sqrt(bs);
				d__1 = -b / a;
				bvn -= exp(-hk / 2) * sqrt(6.283185307179586) * mvphi_(&d__1)
					   * b * (1 - c__ * bs * (1 - d__ * bs / 5) / 3);
			}
			a /= 2;
			i__1 = lg;
			for (i__ = 1; i__ <= i__1; ++i__) {
				/* Computing 2nd power */
				d__1 = a * (x[i__ + ng * 10 - 11] + 1);
				xs = d__1 * d__1;
				rs = sqrt(1 - xs);
				bvn += a * w[i__ + ng * 10 - 11] * (exp(-bs / (xs * 2) - hk /
																		 (rs + 1)) / rs - exp(-(bs / xs + hk) / 2) * (c__ * xs
																													  * (d__ * xs + 1) + 1));
				/* Computing 2nd power */
				d__1 = -x[i__ + ng * 10 - 11] + 1;
				xs = as * (d__1 * d__1) / 4;
				rs = sqrt(1 - xs);
				/* Computing 2nd power */
				d__1 = rs + 1;
				bvn += a * w[i__ + ng * 10 - 11] * exp(-(bs / xs + hk) / 2) *
					   (exp(-hk * xs / (d__1 * d__1 * 2)) / rs - (c__ * xs *
																  (d__ * xs + 1) + 1));
			}
			bvn = -bvn / 6.283185307179586;
		}
		if (*r__ > 0.) {
			d__1 = -max(h__,k);
			bvn += mvphi_(&d__1);
		} else {
			bvn = -bvn;
			if (k > h__) {
				if (h__ < 0.) {
					bvn = bvn + mvphi_(&k) - mvphi_(&h__);
				} else {
					d__1 = -h__;
					d__2 = -k;
					bvn = bvn + mvphi_(&d__1) - mvphi_(&d__2);
				}
			}
		}
	}
	ret_val = bvn;
	return ret_val;
}

#undef x
#undef w


/* Student t Distribution Function */
/*       T */
/* TSTDNT = C   I  ( 1 + y*y/NU )**( -(NU+1)/2 ) dy */
/*   NU -INF */
double mvstdt_(int *nu, double *t)
{
	/* System generated locals */
	double ret_val;

	/* Builtin functions */
	double atan(double), sqrt(double);

	/* Local variables */
	static int j;
	static double rn, ts, tt, csthe, snthe;
	extern double mvphi_(double *);
	static double polyn;

	if (*nu < 1) {
		ret_val = mvphi_(t);
	} else if (*nu == 1) {
		ret_val = (atan(*t) * 2 / 3.141592653589793 + 1) / 2;
	} else if (*nu == 2) {
		ret_val = (*t / sqrt(*t * *t + 2) + 1) / 2;
	} else {
		tt = *t * *t;
		csthe = *nu / (*nu + tt);
		polyn = 1.;
		for (j = *nu - 2; j >= 2; j += -2) {
			polyn = (j - 1) * csthe * polyn / j + 1;
		}
		if (*nu % 2 == 1) {
			rn = (double) (*nu);
			ts = *t / sqrt(rn);
			ret_val = ((atan(ts) + ts * csthe * polyn) * 2 /
					   3.141592653589793 + 1) / 2;
		} else {
			snthe = *t / sqrt(*nu + tt);
			ret_val = (snthe * polyn + 1) / 2;
		}
		if (ret_val < 0.) {
			ret_val = 0.;
		}
	}
	return ret_val;
}

/* A function for computing bivariate normal and t probabilities. */

/* Parameters */

/* NU     INTEGER degrees of freedom parameter; NU < 1 gives normal case. */
/* LOWER  REAL, array of lower integration limits. */
/* UPPER  REAL, array of upper integration limits. */
/* INFIN  INTEGER, array of integration limits ints: */
/* if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; */
/* if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); */
/* if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. */
/* CORREL REAL, correlation coefficient. */
double mvbvt_(int *nu, double *lower, double *upper, int * infin,
			  double *correl)
{
	/* System generated locals */
	double ret_val, d__1, d__2, d__3, d__4;

	/* Local variables */
	extern double mvbvn_(double *, double *, int *,
						 double *), mvbvtl_(int *, double *, double *,
											double *);

	/* Parameter adjustments */
	--infin;
	--upper;
	--lower;

	/* Function Body */
	if (*nu < 1) {
		ret_val = mvbvn_(&lower[1], &upper[1], &infin[1], correl);
	} else {
		if (infin[1] == 2 && infin[2] == 2) {
			ret_val = mvbvtl_(nu, &upper[1], &upper[2], correl) - mvbvtl_(nu,
																		  &upper[1], &lower[2], correl) - mvbvtl_(nu, &lower[1], &
					upper[2], correl) + mvbvtl_(nu, &lower[1], &lower[2],
												correl);
		} else if (infin[1] == 2 && infin[2] == 1) {
			d__1 = -lower[1];
			d__2 = -lower[2];
			d__3 = -upper[1];
			d__4 = -lower[2];
			ret_val = mvbvtl_(nu, &d__1, &d__2, correl) - mvbvtl_(nu, &d__3, &
					d__4, correl);
		} else if (infin[1] == 1 && infin[2] == 2) {
			d__1 = -lower[1];
			d__2 = -lower[2];
			d__3 = -lower[1];
			d__4 = -upper[2];
			ret_val = mvbvtl_(nu, &d__1, &d__2, correl) - mvbvtl_(nu, &d__3, &
					d__4, correl);
		} else if (infin[1] == 2 && infin[2] == 0) {
			ret_val = mvbvtl_(nu, &upper[1], &upper[2], correl) - mvbvtl_(nu,
																		  &lower[1], &upper[2], correl);
		} else if (infin[1] == 0 && infin[2] == 2) {
			ret_val = mvbvtl_(nu, &upper[1], &upper[2], correl) - mvbvtl_(nu,
																		  &upper[1], &lower[2], correl);
		} else if (infin[1] == 1 && infin[2] == 0) {
			d__1 = -lower[1];
			d__2 = -(*correl);
			ret_val = mvbvtl_(nu, &d__1, &upper[2], &d__2);
		} else if (infin[1] == 0 && infin[2] == 1) {
			d__1 = -lower[2];
			d__2 = -(*correl);
			ret_val = mvbvtl_(nu, &upper[1], &d__1, &d__2);
		} else if (infin[1] == 1 && infin[2] == 1) {
			d__1 = -lower[1];
			d__2 = -lower[2];
			ret_val = mvbvtl_(nu, &d__1, &d__2, correl);
		} else if (infin[1] == 0 && infin[2] == 0) {
			ret_val = mvbvtl_(nu, &upper[1], &upper[2], correl);
		} else {
			ret_val = 1.;
		}
	}
	return ret_val;
}

/* A function for computing complementary bivariate normal and t */
/* probabilities. */

/* Parameters */

/* NU     INTEGER degrees of freedom parameter. */
/* L      REAL, array of lower integration limits. */
/* U      REAL, array of upper integration limits. */
/* INFIN  INTEGER, array of integration limits ints: */
/* if INFIN(1) INFIN(2),        then MVBVTC computes */
/* 0         0              P( X>U(1), Y>U(2) ) */
/* 1         0              P( X<L(1), Y>U(2) ) */
/* 0         1              P( X>U(1), Y<L(2) ) */
/* 1         1              P( X<L(1), Y<L(2) ) */
/* 2         0      P( X>U(1), Y>U(2) ) + P( X<L(1), Y>U(2) ) */
/* 2         1      P( X>U(1), Y<L(2) ) + P( X<L(1), Y<L(2) ) */
/* 0         2      P( X>U(1), Y>U(2) ) + P( X>U(1), Y<L(2) ) */
/* 1         2      P( X<L(1), Y>U(2) ) + P( X<L(1), Y<L(2) ) */
/* 2         2      P( X>U(1), Y<L(2) ) + P( X<L(1), Y<L(2) ) */
/*               +  P( X>U(1), Y>U(2) ) + P( X<L(1), Y>U(2) ) */

/* RHO    REAL, correlation coefficient. */
double mvbvtc_(int *nu, double *l, double *u, int *infin,
			   double *rho)
{
	/* System generated locals */
	double ret_val;

	/* Local variables */
	static double b;
	static int i__;
	static double lw[2], up[2];
	static int inf[2];
	extern double mvbvt_(int *, double *, double *, int *,
						 double *);


	/* Parameter adjustments */
	--infin;
	--u;
	--l;

	/* Function Body */
	for (i__ = 1; i__ <= 2; ++i__) {
		if (infin[i__] % 2 == 0) {
			inf[i__ - 1] = 1;
			lw[i__ - 1] = u[i__];
		} else {
			inf[i__ - 1] = 0;
			up[i__ - 1] = l[i__];
		}
	}
	b = mvbvt_(nu, lw, up, inf, rho);
	for (i__ = 1; i__ <= 2; ++i__) {
		if (infin[i__] == 2) {
			inf[i__ - 1] = 0;
			up[i__ - 1] = l[i__];
			b += mvbvt_(nu, lw, up, inf, rho);
		}
	}
	if (infin[1] == 2 && infin[2] == 2) {
		inf[0] = 1;
		lw[0] = u[1];
		b += mvbvt_(nu, lw, up, inf, rho);
	}
	ret_val = b;
	return ret_val;
}

/* a function for computing bivariate t probabilities. */

/* Alan Genz */
/* Department of Mathematics */
/* Washington State University */
/* Pullman, Wa 99164-3113 */
/* Email : alangenz@wsu.edu */

/* this function is based on the method described by */
/* Dunnett, C.W. and M. Sobel, (1954), */
/* A bivariate generalization of Student's t-distribution */
/* with tables for certain special cases, */
/* Biometrika 41, pp. 153-169. */

/* mvbvtl - calculate the probability that x < dh and y < dk. */

/* parameters */

/* nu number of degrees of freedom */
/* dh 1st lower integration limit */
/* dk 2nd lower integration limit */
/* r   correlation coefficient */
double mvbvtl_(int *nu, double *dh, double *dk, double *r__)
{
	/* System generated locals */
	int i__1;
	double ret_val, d__1, d__2, d__3;

	/* Builtin functions */
	double sqrt(double), atan2(double, double);

	/* Local variables */
	static int j, hs, ks;
	static double hkn, hpk, hrk, krh, bvt, ors, snu, gmph, gmpk, hkrn,
			qhrk, xnkh, xnhk, btnckh, btnchk, btpdkh, btpdhk;

	snu = sqrt((double) (*nu));
	ors = 1 - *r__ * *r__;
	hrk = *dh - *r__ * *dk;
	krh = *dk - *r__ * *dh;
	if (fabs(hrk) + ors > 0.) {
		/* Computing 2nd power */
		d__1 = hrk;
		/* Computing 2nd power */
		d__2 = hrk;
		/* Computing 2nd power */
		d__3 = *dk;
		xnhk = d__1 * d__1 / (d__2 * d__2 + ors * (*nu + d__3 * d__3));
		/* Computing 2nd power */
		d__1 = krh;
		/* Computing 2nd power */
		d__2 = krh;
		/* Computing 2nd power */
		d__3 = *dh;
		xnkh = d__1 * d__1 / (d__2 * d__2 + ors * (*nu + d__3 * d__3));
	} else {
		xnhk = 0.;
		xnkh = 0.;
	}
	d__1 = *dh - *r__ * *dk;
	hs = (int) d_sign(&c_b91, &d__1);
	d__1 = *dk - *r__ * *dh;
	ks = (int) d_sign(&c_b91, &d__1);
	if (*nu % 2 == 0) {
		bvt = atan2(sqrt(ors), -(*r__)) / 6.2831853071795862;
		/* Computing 2nd power */
		d__1 = *dh;
		gmph = *dh / sqrt((*nu + d__1 * d__1) * 16);
		/* Computing 2nd power */
		d__1 = *dk;
		gmpk = *dk / sqrt((*nu + d__1 * d__1) * 16);
		btnckh = atan2(sqrt(xnkh), sqrt(1 - xnkh)) * 2 /
				 3.14159265358979323844;
		btpdkh = sqrt(xnkh * (1 - xnkh)) * 2 / 3.14159265358979323844;
		btnchk = atan2(sqrt(xnhk), sqrt(1 - xnhk)) * 2 /
				 3.14159265358979323844;
		btpdhk = sqrt(xnhk * (1 - xnhk)) * 2 / 3.14159265358979323844;
		i__1 = *nu / 2;
		for (j = 1; j <= i__1; ++j) {
			bvt += gmph * (ks * btnckh + 1);
			bvt += gmpk * (hs * btnchk + 1);
			btnckh += btpdkh;
			btpdkh = (j << 1) * btpdkh * (1 - xnkh) / ((j << 1) + 1);
			btnchk += btpdhk;
			btpdhk = (j << 1) * btpdhk * (1 - xnhk) / ((j << 1) + 1);
			/* Computing 2nd power */
			d__1 = *dh;
			gmph = gmph * ((j << 1) - 1) / ((j << 1) * (d__1 * d__1 / *nu + 1)
			);
			/* Computing 2nd power */
			d__1 = *dk;
			gmpk = gmpk * ((j << 1) - 1) / ((j << 1) * (d__1 * d__1 / *nu + 1)
			);
		}
	} else {
		/* Computing 2nd power */
		d__1 = *dh;
		/* Computing 2nd power */
		d__2 = *dk;
		qhrk = sqrt(d__1 * d__1 + d__2 * d__2 - *r__ * 2 * *dh * *dk + *nu * ors);
		hkrn = *dh * *dk + *r__ * *nu;
		hkn = *dh * *dk - *nu;
		hpk = *dh + *dk;
		bvt = atan2(-snu * (hkn * qhrk + hpk * hkrn), hkn * hkrn - *nu * hpk *
																   qhrk) / 6.2831853071795862;
		if (bvt < -1e-15) {
			bvt += 1;
		}
		/* Computing 2nd power */
		d__1 = *dh;
		gmph = *dh / (snu * 6.2831853071795862 * (d__1 * d__1 / *nu + 1));
		/* Computing 2nd power */
		d__1 = *dk;
		gmpk = *dk / (snu * 6.2831853071795862 * (d__1 * d__1 / *nu + 1));
		btnckh = sqrt(xnkh);
		btpdkh = btnckh;
		btnchk = sqrt(xnhk);
		btpdhk = btnchk;
		i__1 = (*nu - 1) / 2;
		for (j = 1; j <= i__1; ++j) {
			bvt += gmph * (ks * btnckh + 1);
			bvt += gmpk * (hs * btnchk + 1);
			btpdkh = ((j << 1) - 1) * btpdkh * (1 - xnkh) / (j << 1);
			btnckh += btpdkh;
			btpdhk = ((j << 1) - 1) * btpdhk * (1 - xnhk) / (j << 1);
			btnchk += btpdhk;
			/* Computing 2nd power */
			d__1 = *dh;
			gmph = (j << 1) * gmph / (((j << 1) + 1) * (d__1 * d__1 / *nu + 1)
			);
			/* Computing 2nd power */
			d__1 = *dk;
			gmpk = (j << 1) * gmpk / (((j << 1) + 1) * (d__1 * d__1 / *nu + 1)
			);
		}
	}
	ret_val = bvt;
	return ret_val;
}


/*  MVCHNV */
/* P =  1 - K  I     exp(-t*t/2) t**(N-1) dt, for N >= 1. */
/* N  0 */
double mvchnv_(int *n, double *p)
{
	/* Initialized data */

	static int no = 0;

	/* System generated locals */
	double ret_val, d__1;

	/* Builtin functions */
	double log(double), sqrt(double), exp(double);

	/* Local variables */
	static int i__;
	static double r__, ro, lkn;
	extern double mvchnc_(double *, int *, double *,
						  double *), mvphnv_(double *);

	/* LRP =   LOG( SQRT( 2/PI ) ) */
	if (*n <= 1) {
		d__1 = *p / 2;
		r__ = -mvphnv_(&d__1);
	} else if (*p < 1.) {
		if (*n == 2) {
			r__ = sqrt(log(*p) * -2);
		} else {
			if (*n != no) {
				no = *n;
				lkn = 0.;
				for (i__ = *n - 2; i__ >= 2; i__ += -2) {
					lkn -= log((double) i__);
				}
				if (*n % 2 == 1) {
					lkn += -.22579135264472743235;
				}
			}
			if ((double) (*n) >= log(1 - *p) * -5 / 4) {
				r__ = 2. / (*n * 9);
				/* Computing 3rd power */
				d__1 = -mvphnv_(p) * sqrt(r__) + 1 - r__;
				r__ = *n * (d__1 * (d__1 * d__1));
				if (r__ > (double) ((*n << 1) + 6)) {
					r__ = (lkn - log(*p)) * 2 + (*n - 2) * log(r__);
				}
			} else {
				r__ = exp((log((1 - *p) * *n) - lkn) * 2. / *n);
			}
			r__ = sqrt(r__);
			ro = r__;
			r__ = mvchnc_(&lkn, n, p, &r__);
			if ((d__1 = r__ - ro, fabs(d__1)) > 1e-6) {
				ro = r__;
				r__ = mvchnc_(&lkn, n, p, &r__);
				if ((d__1 = r__ - ro, fabs(d__1)) > 1e-6) {
					r__ = mvchnc_(&lkn, n, p, &r__);
				}
			}
		}
	} else {
		r__ = 0.;
	}
	ret_val = r__;
	return ret_val;
}

/* Third order Schroeder correction to R for MVCHNV */
double mvchnc_(double *lkn, int *n, double *p, double *r__)
{
	/* System generated locals */
	double ret_val, d__1;

	/* Builtin functions */
	double log(double), exp(double);

	/* Local variables */
	static int i__;
	static double df, ai, bi, al, ci, di, dl, rn, rr, chi;
	extern double mvphi_(double *);

	/* LRP =   LOG( SQRT( 2/PI ) ) */
	rr = *r__ * *r__;
	if (*n < 2) {
		d__1 = -(*r__);
		chi = mvphi_(&d__1) * 2;
	} else if (*n < 100) {

		/* Use standard Chi series */
		rn = 1.;
		for (i__ = *n - 2; i__ >= 2; i__ += -2) {
			rn = rr * rn / i__ + 1;
		}
		rr /= 2;
		if (*n % 2 == 0) {
			chi = exp(log(rn) - rr);
		} else {
			d__1 = -(*r__);
			chi = exp(log(*r__ * rn) - .22579135264472743235 - rr) + mvphi_(&
																					d__1) * 2;
		}
	} else {
		rr /= 2;
		al = *n / 2.;
		chi = exp(-rr + al * log(rr) + *lkn + log(2.) * (*n - 2) / 2);
		if (rr < al + 1) {

			/* Use Incomplete Gamma series */
			dl = chi;
			for (i__ = 1; i__ <= 1000; ++i__) {
				dl = dl * rr / (al + i__);
				chi += dl;
				if ((d__1 = dl * rr / (al + i__ + 1 - rr), fabs(d__1)) < 1e-14)
				{
					goto L10;
				}
			}
			L10:
			chi = 1 - chi / al;
		} else {

			/* Use Incomplete Gamma continued fraction */
			bi = rr + 1 - al;
			ci = 1e14;
			di = bi;
			chi /= bi;
			for (i__ = 1; i__ <= 250; ++i__) {
				ai = i__ * (al - i__);
				bi += 2;
				ci = bi + ai / ci;
				if (ci == 0.) {
					ci = 1e-14;
				}
				di = bi + ai / di;
				if (di == 0.) {
					di = 1e-14;
				}
				dl = ci / di;
				chi *= dl;
				if ((d__1 = dl - 1, fabs(d__1)) < 1e-14) {
					goto L20;
				}
			}
		}
	}
	L20:
	df = (*p - chi) / exp(*lkn + (*n - 1) * log(*r__) - rr);
	ret_val = *r__ - df * (1 - df * (*r__ - (*n - 1) / *r__) / 2);
	return ret_val;
}

/* Automatic Multidimensional Integration Subroutine */

/* AUTHOR: Alan Genz */
/* Department of Mathematics */
/* Washington State University */
/* Pulman, WA 99164-3113 */
/* Email: AlanGenz@wsu.edu */

/* Last Change: 12/15/00 */

/* MVKBRV computes an approximation to the integral */

/* 1  1     1 */
/* I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1) */
/* 0  0     0 */

/* F(X) is a float NF-vector of integrands. */

/* It uses randomized Korobov rules. The primary references are */
/* "Randomization of Number Theoretic Methods for Multiple Integration" */
/* R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14, */
/* and */
/* "Optimal Parameters for Multidimensional Integration", */
/* P. Keast, SIAM J Numer Anal, 10, pp.831-838. */
/* If there are more than 100 variables, the remaining variables are */
/* integrated using the rules described in the reference */
/* "On a Number-Theoretical Integration Method" */
/* H. Niederreiter, Aequationes Mathematicae, 8(1972), pp. 304-11. */

/* **************  Parameters ******************************************** */
/* ***** Input parameters */
/* NDIM    Number of variables, must exceed 1, but not exceed 100 */
/* MINVLS  Integer minimum number of function evaluations allowed. */
/* MINVLS must not exceed MAXVLS.  If MINVLS < 0 then the */
/* routine assumes a previous call has been made with */
/* the same integrands and continues that calculation. */
/* MAXVLS  Integer maximum number of function evaluations allowed. */
/* NF      Number of integrands, must exceed 1, but not exceed 5000 */
/* FUNSUB  EXTERNALly declared user defined integrand subroutine. */
/* It must have parameters ( NDIM, Z, NF, FUNVLS ), where */
/* Z is a float NDIM-vector and FUNVLS is a float NF-vector. */

/* ABSEPS  Required absolute accuracy. */
/* RELEPS  Required relative accuracy. */
/* ***** Output parameters */
/* MINVLS  Actual number of function evaluations used. */
/* ABSERR  Maximum norm of estimated absolute accuracy of FINEST. */
/* FINEST  Estimated NF-vector of values of the integrals. */
/* INFORM  INFORM = 0 for normal exit, when */
/*     ABSERR <= MAX(ABSEPS, RELEPS*||FINEST||) */
/*  and */
/*     INTVLS <= MAXCLS. */
/* INFORM = 1 If MAXVLS was too small to obtain the required */
/* accuracy. In this case a value FINEST is returned with */
/* estimated absolute accuracy ABSERR. */
/* *********************************************************************** */
int mvkbrv_(int *ndim, int *minvls, int *maxvls,
			int *nf, U_fp funsub, double *abseps, double *releps,
			double *abserr, double *finest, int *inform__)
{
	/* Initialized data */
	static int p[28] = { 31,47,73,113,173,263,397,593,907,1361,2053,3079,
							  4621,6947,10427,15641,23473,35221,52837,79259,118891,178349,
							  267523,401287,601943,902933,1354471,2031713 };
	static int c__[2772]	/* was [28][99] */ = { 12,13,27,35,64,111,163,
														246,347,505,794,1189,1763,2872,4309,6610,9861,10327,19540,34566,
														31929,40701,103650,165843,130365,333459,500884,858339,9,11,28,27,
														66,42,154,189,402,220,325,888,1018,3233,3758,6977,3647,7582,19926,
														9579,49367,69087,125480,90647,236711,375354,566009,918142,9,17,10,
														27,28,54,83,242,322,601,960,259,1500,1534,4034,1686,4073,7124,
														11582,12654,10982,77576,59978,59925,110235,102417,399251,501970,
														13,10,11,36,28,118,43,102,418,644,528,1082,432,2941,1963,3819,
														2535,8214,11113,26856,3527,64590,46875,189541,125699,383544,
														652979,234813,12,15,11,22,44,20,82,250,215,612,247,725,1332,2910,
														730,2314,3430,9600,24585,37873,27066,39397,77172,67647,56483,
														292630,355008,460565,12,15,20,29,44,31,92,250,220,160,247,811,
														2203,393,642,5647,9865,10271,8726,38806,13226,33179,83021,74795,
														93735,41147,430235,31996,12,15,11,29,55,31,150,102,339,206,338,
														636,126,1796,1502,3953,2830,10193,17218,29501,56010,10858,126904,
														68365,234469,374614,328722,753018,12,15,11,20,67,72,59,250,339,
														206,366,965,2240,919,2246,3614,9328,10800,419,17271,18911,38935,
														14541,167485,60549,48032,670680,256150,12,15,28,45,10,17,76,280,
														339,206,847,497,1719,446,3834,5115,4320,9086,4918,3663,40574,
														43129,56299,143918,1291,435453,405585,199809,12,15,13,5,10,94,76,
														118,337,422,753,497,1284,919,1511,423,5913,2365,4918,10763,20767,
														35468,43636,74912,93937,281493,405585,993599,12,22,13,5,10,14,47,
														196,218,134,753,1490,878,919,1102,423,10365,4409,4918,18955,20767,
														35468,11655,167289,245291,358168,424646,245149,12,15,28,5,10,14,
														11,118,315,518,236,1490,1983,1117,1102,5408,8272,13812,15701,1298,
														9686,5279,52680,75517,196061,114121,670180,794183,3,15,13,21,10,
														11,11,191,315,134,334,392,266,103,1522,7426,3706,5661,17710,26560,
														47603,61518,88549,8148,258647,346892,670180,121349,3,6,13,21,10,
														14,100,215,315,134,334,1291,266,103,1522,423,6186,9344,4037,17132,
														47603,61518,29804,172106,162489,238990,641587,150619,3,6,13,21,38,
														14,131,121,315,518,461,508,266,103,3427,423,7806,9344,4037,17132,
														11736,27945,101894,126159,176631,317313,215580,376952,12,6,14,21,
														38,14,116,121,167,652,711,508,266,103,3427,487,7806,10362,15808,
														4753,11736,70975,113675,35867,204895,164158,59048,809123,7,15,14,
														21,10,94,116,49,167,382,652,1291,747,103,3928,6227,7806,9344,
														11401,4753,41601,70975,48040,35867,73353,35497,633320,809123,7,15,
														14,21,10,10,116,49,167,206,381,1291,747,103,915,2660,8610,9344,
														19398,8713,12888,86478,113675,35867,172319,70530,81010,804319,12,
														9,14,21,10,10,116,49,167,158,381,508,127,103,915,6227,2563,8585,
														25950,18624,32948,86478,34987,121694,28881,70530,20789,67352,12,
														13,14,21,10,10,116,49,361,441,381,1291,127,2311,3818,1221,11558,
														11114,25950,13082,30801,20514,48308,52171,136787,434839,389250,
														969594,12,2,14,21,10,10,116,49,201,179,652,508,2074,3117,3818,
														3811,11558,13080,4454,6791,44243,20514,97926,95354,122081,24754,
														389250,434796,12,2,14,21,49,14,138,49,124,441,381,508,127,1101,
														3818,197,9421,13080,24987,1122,53351,73178,5475,113969,122081,
														24754,638764,969594,12,2,14,21,49,14,138,49,124,56,381,867,2074,
														3117,3818,4367,1181,13080,11719,19363,53351,73178,49449,113969,
														275993,24754,638764,804319,12,13,14,21,49,14,138,49,124,559,381,
														867,1400,3117,4782,351,9421,6949,8697,34695,16016,43098,6850,
														76304,64673,393656,389250,391368,12,11,14,21,49,14,138,49,124,559,
														381,867,1383,1101,4782,1281,1181,3436,1452,18770,35086,43098,
														62545,123709,211587,118711,389250,761041,12,11,14,21,49,14,138,49,
														124,56,381,867,1383,1101,4782,1221,1181,3436,1452,18770,35086,
														4701,62545,123709,211587,118711,398094,754049,12,10,14,21,49,14,
														138,49,124,56,381,934,1383,1101,3818,351,1181,3436,1452,18770,
														32581,59979,9440,144615,211587,148227,80846,466264,3,15,14,21,49,
														14,138,49,124,56,381,867,1383,1101,4782,351,9421,13213,1452,18770,
														2464,59979,33242,123709,282859,271087,147776,754049,3,15,14,29,49,
														11,138,171,124,56,226,867,1383,1101,3818,351,1181,6130,1452,15628,
														2464,58556,9440,64958,282859,355831,147776,754049,3,15,14,17,49,
														11,138,171,124,56,326,867,1383,2503,3818,7245,1181,6130,8697,
														18770,49554,69916,33242,64958,211587,91034,296177,466264,12,15,14,
														17,49,11,101,171,124,56,326,867,1383,2503,1327,1984,10574,8159,
														8697,18770,2464,15170,9440,32377,242821,417029,398094,754049,7,15,
														31,17,49,8,101,171,124,56,326,867,1383,2503,1327,2999,10574,8159,
														6436,18770,2464,15170,33242,193002,256865,417029,398094,754049,7,
														15,31,17,49,8,101,171,231,56,326,867,1383,2503,1327,2999,3534,
														11595,21475,18770,49554,4832,9440,193002,256865,91034,147776,
														282852,12,15,5,17,38,8,101,171,231,56,326,867,1383,2503,1327,2999,
														3534,8159,6436,33766,49554,4832,62850,25023,256865,91034,147776,
														429907,12,15,5,17,38,8,101,171,90,56,326,1284,1400,2503,1327,2999,
														3534,3436,22913,20837,2464,43064,9440,40017,122203,417029,396313,
														390017,12,15,5,17,31,8,101,171,90,56,326,1284,1383,2503,1327,2999,
														3534,7096,6434,20837,81,71685,9440,141605,291915,91034,578233,
														276645,12,6,31,17,4,8,101,171,90,56,126,1284,1383,2503,1327,2999,
														3534,7096,18497,20837,27260,4832,9440,189165,122203,299843,578233,
														994856,12,6,13,17,4,8,101,171,90,56,326,1284,1383,429,1387,3995,
														2898,7096,11089,20837,10681,15170,90308,189165,291915,299843,
														578233,250142,12,6,11,17,31,18,101,171,90,56,326,1284,1383,429,
														1387,2063,2898,7096,11089,20837,2185,15170,90308,141605,291915,
														413548,19482,144595,12,15,11,23,64,18,101,171,90,101,326,1284,
														1383,429,1387,2063,2898,7096,11089,20837,2185,15170,90308,189165,
														122203,413548,620706,907454,12,15,11,23,4,18,101,171,90,101,326,
														1284,1383,429,1387,2063,3450,7096,11089,6545,2185,27679,47904,
														189165,25639,308300,187095,689648,12,9,11,23,4,18,101,171,90,56,
														326,1284,1383,429,1387,2063,2141,7096,3036,6545,2185,27679,47904,
														141605,25639,413548,620706,687580,3,13,11,23,4,18,101,171,90,101,
														326,1284,507,429,1387,1644,2141,7096,3036,6545,2185,27679,47904,
														141605,291803,413548,187095,687580,3,2,11,23,64,113,101,171,90,
														101,326,563,1073,429,1387,2063,2141,7096,14208,6545,2185,60826,
														47904,141605,245397,413548,126467,687580,3,2,13,23,45,62,101,171,
														90,101,326,563,1073,1702,1387,2077,2141,7096,14208,6545,2185,
														60826,47904,189165,284047,308300,241663,687580,12,2,13,23,45,62,
														101,171,90,101,326,563,1073,1702,1387,2512,2141,7096,14208,12138,
														18086,6187,47904,127047,245397,308300,241663,978368,7,13,13,23,45,
														45,101,171,90,101,326,563,1073,1702,2339,2512,2141,7096,14208,
														12138,18086,6187,47904,127047,245397,308300,241663,687580,7,11,13,
														23,45,45,101,171,90,101,195,1010,1990,184,2339,2512,2141,7096,
														12906,12138,18086,4264,47904,127047,245397,413548,241663,552742,
														12,11,13,23,45,113,101,171,48,101,195,1010,1990,184,2339,2077,
														7055,7096,12906,12138,18086,4264,47904,127047,245397,308300,
														241663,105195,12,10,13,23,45,113,101,171,48,101,55,1010,1990,184,
														2339,2077,7055,7096,12906,12138,18086,4264,41143,127047,245397,
														308300,241663,942843,12,15,13,23,66,113,101,171,48,193,55,208,
														1990,184,2339,2077,7055,7096,12906,12138,17631,4264,41143,127047,
														245397,308300,241663,768249,12,15,14,21,66,113,116,171,48,193,55,
														838,1990,184,2339,2077,7055,7096,12906,12138,17631,4264,41143,
														127047,245397,308300,241663,307142,12,15,14,27,66,113,116,171,90,
														193,55,563,507,105,2339,754,7055,7096,12906,12138,18086,45567,
														41143,127047,94241,308300,241663,307142,12,15,14,3,66,113,116,171,
														90,193,55,563,507,105,2339,754,7055,4377,12906,12138,18086,32269,
														41143,127047,66575,15311,241663,307142,12,15,14,3,66,113,116,171,
														90,193,55,563,507,105,2339,754,7055,7096,12906,12138,18086,32269,
														41143,127047,66575,15311,241663,307142,12,15,14,3,66,113,116,171,
														90,193,55,759,507,105,2339,754,7055,4377,7614,12138,37335,32269,
														41143,127047,217673,15311,241663,880619,12,15,14,24,66,113,116,
														171,90,193,55,759,507,105,2339,754,7055,4377,7614,12138,37774,
														32269,36114,127047,217673,15311,321632,880619,3,15,14,27,66,113,
														100,171,90,101,55,564,507,105,2339,754,7055,4377,7614,12138,37774,
														62060,36114,127047,217673,176255,23210,880619,3,15,14,27,66,113,
														100,171,90,101,55,759,507,105,2339,754,7055,4377,7614,12138,37774,
														62060,36114,127047,217673,176255,23210,880619,3,6,14,17,66,113,
														100,171,90,101,55,759,507,105,3148,754,7055,4377,5021,30483,26401,
														62060,36114,127047,217673,23613,394484,880619,12,6,14,29,66,113,
														100,171,90,101,55,801,507,105,3148,754,7055,5410,5021,30483,26401,
														62060,36114,127047,217673,23613,394484,880619,7,6,14,29,66,113,
														100,171,90,101,55,801,1073,105,3148,754,7055,5410,5021,30483,
														26401,62060,24997,127047,217673,23613,394484,880619,7,15,14,29,66,
														113,138,161,90,101,55,801,1073,105,3148,754,7055,4377,5021,30483,
														26401,62060,65162,127047,217673,23613,78101,117185,12,15,14,17,66,
														113,138,161,90,101,55,801,1073,105,3148,754,2831,4377,5021,30483,
														26401,62060,65162,127047,217673,23613,78101,117185,12,9,14,5,66,
														113,138,161,90,101,55,759,1073,105,3148,754,8204,4377,5021,12138,
														26401,62060,65162,127047,217673,23613,78101,117185,12,13,14,5,66,
														63,138,161,90,101,55,759,1073,105,3148,754,8204,4377,10145,12138,
														26401,62060,65162,127785,217673,172210,542095,117185,12,2,14,5,66,
														63,138,161,90,101,55,759,1073,105,3148,754,8204,4377,10145,12138,
														26401,1803,65162,127785,217673,204328,542095,117185,12,2,31,5,66,
														53,101,161,90,101,55,759,1073,105,3148,754,8204,4377,10145,12138,
														26401,1803,65162,127785,217673,204328,542095,117185,12,2,31,21,66,
														63,101,161,90,101,195,759,1073,105,3148,754,8204,4377,10145,12138,
														26401,1803,65162,127785,217673,204328,542095,117185,12,13,5,21,11,
														67,101,161,90,101,195,563,1073,105,3148,754,8204,4377,10145,12138,
														26401,1803,65162,127785,217673,204328,542095,117185,12,11,5,21,66,
														67,101,14,90,101,195,563,1073,105,3148,754,8204,4377,10145,12138,
														26401,1803,65162,127785,217673,121626,542095,117185,7,11,5,21,66,
														67,101,14,90,101,195,563,1073,105,3148,1097,8204,4377,10145,12138,
														26401,1803,65162,127785,217673,121626,542095,117185,3,10,11,21,66,
														67,101,14,90,101,195,563,1073,105,3148,1097,8204,4377,10145,12138,
														12982,1803,65162,127785,217673,121626,542095,117185,3,10,13,21,66,
														67,101,14,90,101,195,563,1073,105,3148,754,8204,4377,10145,12138,
														40398,1803,65162,127785,217673,121626,542095,60731,3,15,11,21,66,
														67,101,14,90,101,195,563,1073,105,3148,754,8204,4377,10145,12138,
														40398,1803,65162,127785,210249,121626,542095,60731,7,15,11,21,66,
														67,101,14,243,101,132,563,1073,105,3148,754,8204,4377,10145,12138,
														40398,1803,65162,80822,210249,200187,542095,60731,7,15,11,21,66,
														67,101,14,243,101,132,563,1073,105,3148,754,8204,4377,10145,12138,
														40398,1803,47650,80822,210249,200187,542095,60731,7,15,11,21,66,
														67,101,14,243,101,132,226,1073,105,1776,248,8204,4377,10145,12138,
														40398,1803,47650,80822,210249,200187,542095,60731,3,15,11,21,66,
														67,101,14,243,122,132,226,22,105,1776,754,8204,4377,10145,12138,
														40398,1803,47650,80822,210249,200187,542095,60731,3,15,11,21,45,
														67,101,14,243,122,132,226,22,105,1776,1097,8204,4377,10145,12138,
														3518,51108,47650,80822,210249,200187,542095,60731,3,15,11,21,11,
														67,101,14,243,122,132,226,22,105,3354,1097,8204,4377,10145,12138,
														3518,51108,47650,80822,210249,121551,542095,60731,3,15,13,21,7,67,
														101,14,243,122,132,226,22,105,3354,1097,8204,4377,10145,12138,
														3518,51108,47650,131661,210249,121551,542095,60731,3,6,13,21,3,67,
														101,14,243,122,132,226,22,105,3354,1097,8204,4377,10145,12138,
														37799,51108,47650,131661,210249,248492,542095,60731,3,2,11,21,2,
														67,101,14,243,122,132,226,22,105,925,222,8204,4377,10145,9305,
														37799,51108,40586,131661,210249,248492,542095,60731,3,3,13,17,2,
														51,101,14,243,122,132,226,1073,105,3354,222,8204,4377,10145,11107,
														37799,51108,40586,131661,94453,248492,277743,178309,3,2,5,17,2,51,
														101,14,283,122,132,226,452,105,3354,222,8204,4377,10145,11107,
														37799,51108,40586,131661,94453,248492,277743,178309,3,3,5,17,27,
														51,38,14,283,122,387,226,452,784,925,222,8204,4377,10145,11107,
														37799,51108,40586,131661,94453,248492,277743,178309,3,2,5,6,5,51,
														38,10,283,122,387,226,452,784,925,754,8204,4377,10145,11107,37799,
														51108,40586,131661,94453,248492,457259,178309,3,2,5,17,3,51,38,10,
														283,122,387,226,452,784,925,1982,4688,4377,10145,11107,37799,
														51108,40586,131661,94453,248492,457259,74373,3,2,14,17,3,12,38,10,
														283,122,387,226,452,784,925,1982,4688,4377,4544,11107,37799,51108,
														40586,131661,94453,248492,457259,74373,3,2,13,6,5,51,38,10,283,
														122,387,226,452,784,925,1982,4688,4377,4544,11107,37799,51108,
														38725,131661,94453,248492,457259,74373,3,2,5,3,5,12,38,10,283,122,
														387,226,318,784,2133,1982,2831,4377,4544,11107,4721,55315,38725,
														131661,94453,248492,457259,74373,3,2,5,6,2,51,38,10,283,122,387,
														226,301,784,2133,1982,2831,4377,4544,11107,4721,55315,38725,
														131661,94453,248492,457259,74373,3,2,5,6,2,5,38,103,283,122,387,
														226,301,784,2133,1982,2831,4377,4544,11107,4721,54140,38725,
														131661,94453,248492,457259,74373,3,2,5,3,2,3,3,10,16,122,387,226,
														301,784,2133,1982,2831,440,4544,11107,4721,54140,88329,131661,
														94453,13942,457259,74373,3,2,5,3,2,3,3,10,283,101,387,226,301,784,
														2133,1982,2831,440,8394,11107,7067,54140,88329,131661,94453,13942,
														457259,74373,3,2,5,3,2,2,3,10,16,101,387,226,86,784,2133,1982,
														2831,1199,8394,11107,7067,54140,88329,131661,94453,13942,457259,
														214965,3,2,5,3,2,2,3,10,283,101,387,226,86,784,2133,1982,2831,
														1199,8394,9305,7067,54140,88329,7114,94453,13942,457259,214965,3,
														2,5,3,2,5,3,5,283,101,387,226,15,784,2133,1982,2831,1199,8394,
														9305,7067,13134,88329,131661,94453,13942,457259,214965 };

	/* System generated locals */
	int i__1, i__2;
	double d__1, d__2, d__3;

	/* Builtin functions */
	double sqrt(double);

	/* Local variables */
	static int i__, k;
	static double r__[1000], x[1000], fs[5000];
	static int np;
	static double vk[1000];
	static int pr[1000], kmx;
	static double difint, finval[5000], varprd;
	static int sampls;
	static double values[5000], varest[5000], varsqr[5000];
	static int intvls;
	extern int mvkrsv_(int *, int *, double *,
					   int *, double *, int *, U_fp, double *,
					   double *, int *, double *);

	/* Parameter adjustments */
	--finest;

	/* Function Body */
	*inform__ = 1;
	intvls = 0;
	varprd = 0.;
	if (*minvls >= 0) {
		i__1 = *nf;
		for (k = 1; k <= i__1; ++k) {
			finest[k] = 0.;
			varest[k - 1] = 0.;
		}
		sampls = 8;
		for (i__ = min(*ndim,10); i__ <= 28; ++i__) {
			np = i__;
			if (*minvls < (sampls << 1) * p[i__ - 1]) {
				goto L10;
			}
		}
		/* Computing MAX */
		i__1 = 8, i__2 = *minvls / (p[np - 1] << 1);
		sampls = max(i__1,i__2);
	}
	L10:
	vk[0] = 1. / p[np - 1];
	k = 1;
	i__1 = *ndim;
	for (i__ = 2; i__ <= i__1; ++i__) {
		if (i__ <= 100) {
			/* Computing MIN */
			i__2 = *ndim - 1;
			d__1 = c__[np + min(i__2,99) * 28 - 29] * (double) k;
			d__2 = (double) p[np - 1];
			k = (int) fmod(d__1, d__2);
			vk[i__ - 1] = k * vk[0];
		} else {
			d__1 = (double) (i__ - 100) / (*ndim - 99);
			vk[i__ - 1] = (double) ((int) (p[np - 1] * pow(c_b99,d__1)));
			d__1 = vk[i__ - 1] / p[np - 1];
			vk[i__ - 1] = fmod(d__1, c_b91);
		}
	}
	i__1 = *nf;
	for (k = 1; k <= i__1; ++k) {
		finval[k - 1] = 0.;
		varsqr[k - 1] = 0.;
	}

	i__1 = sampls;
	for (i__ = 1; i__ <= i__1; ++i__) {
		mvkrsv_(ndim, &c__100, values, &p[np - 1], vk, nf, (U_fp)funsub, x,
				r__, pr, fs);
		i__2 = *nf;
		for (k = 1; k <= i__2; ++k) {
			difint = (values[k - 1] - finval[k - 1]) / i__;
			finval[k - 1] += difint;
			/* Computing 2nd power */
			d__1 = difint;
			varsqr[k - 1] = (i__ - 2) * varsqr[k - 1] / i__ + d__1 * d__1;
		}
	}

	intvls += (sampls << 1) * p[np - 1];
	kmx = 1;
	i__1 = *nf;
	for (k = 1; k <= i__1; ++k) {
		varprd = varest[k - 1] * varsqr[k - 1];
		finest[k] += (finval[k - 1] - finest[k]) / (varprd + 1);
		if (varsqr[k - 1] > 0.) {
			varest[k - 1] = (varprd + 1) / varsqr[k - 1];
		}
		if ((d__1 = finest[k], fabs(d__1)) > (d__2 = finest[kmx], fabs(d__2))) {
			kmx = k;
		}
	}
	*abserr = sqrt(varsqr[kmx - 1] / (varprd + 1)) * 7 / 2;
	/* Computing MAX */
	d__2 = *abseps, d__3 = (d__1 = finest[kmx], fabs(d__1)) * *releps;
	if (*abserr > max(d__2,d__3)) {
		if (np < 28) {
			++np;
		} else {
			/* Computing MIN */
			i__1 = sampls * 3 / 2, i__2 = (*maxvls - intvls) / (p[np - 1] <<
																		  1);
			sampls = min(i__1,i__2);
			sampls = max(8,sampls);
		}
		if (intvls + (sampls << 1) * p[np - 1] <= *maxvls) {
			goto L10;
		}
	} else {
		*inform__ = 0;
	}
	*minvls = intvls;

	return 0;
}

/* For lattice rule sums */
int mvkrsv_(int *ndim, int *kl, double *values,
			int *prime, double *vk, int *nf, U_fp funsub, double *
x, double *r__, int *pr, double *fs)
{
	/* System generated locals */
	int i__1, i__2;
	double d__1;

	/* Local variables */
	static int j, k, jp;
	extern double mvuni_(void);

	/* Parameter adjustments */
	--fs;
	--pr;
	--r__;
	--x;
	--vk;
	--values;

	/* Function Body */
	i__1 = *nf;
	for (j = 1; j <= i__1; ++j) {
		values[j] = 0.;
	}

	/* Determine random shifts for each variable; scramble lattice rule */
	i__1 = *ndim;
	for (j = 1; j <= i__1; ++j) {
		r__[j] = mvuni_();
		if (j < *kl) {
			jp = (int) (j * r__[j] + 1);
			if (jp < j) {
				pr[j] = pr[jp];
			}
			pr[jp] = j;
		} else {
			pr[j] = j;
		}
	}

	/* Compute latice rule sums */
	i__1 = *prime;
	for (k = 1; k <= i__1; ++k) {
		i__2 = *ndim;
		for (j = 1; j <= i__2; ++j) {
			r__[j] += vk[pr[j]];
			if (r__[j] > 1.) {
				--r__[j];
			}
			x[j] = (d__1 = r__[j] * 2 - 1, fabs(d__1));
		}
		(*funsub)(ndim, &x[1], nf, &fs[1]);
		i__2 = *nf;
		for (j = 1; j <= i__2; ++j) {
			values[j] += (fs[j] - values[j]) / ((k << 1) - 1);
		}
		i__2 = *ndim;
		for (j = 1; j <= i__2; ++j) {
			x[j] = 1 - x[j];
		}
		(*funsub)(ndim, &x[1], nf, &fs[1]);
		i__2 = *nf;
		for (j = 1; j <= i__2; ++j) {
			values[j] += (fs[j] - values[j]) / (k << 1);
		}
	}

	return 0;
}

/* Uniform (0,1) random number generator */

/* Reference: */
/* L'Ecuyer, Pierre (1996), */
/* "Combined Multiple Recursive Random Number Generators" */
/* Operations Research 44, pp. 816-822. */
double mvuni_(void)
{
	/* Initialized data */

	static int x10 = 15485857;
	static int x11 = 17329489;
	static int x12 = 36312197;
	static int x20 = 55911127;
	static int x21 = 75906931;
	static int x22 = 96210113;

	/* System generated locals */
	double ret_val;

	/* Local variables */
	static int h__, z__, p12, p13, p21, p23;

	/* INVMP1 = 1/( M1 + 1 ) */

	/* Some eight digit primes for seeds */


	/* Component 1 */
	h__ = x10 / 11714;
	p13 = (x10 - h__ * 11714) * 183326 - h__ * 2883;
	h__ = x11 / 33921;
	p12 = (x11 - h__ * 33921) * 63308 - h__ * 12979;
	if (p13 < 0) {
		p13 += 2147483647;
	}
	if (p12 < 0) {
		p12 += 2147483647;
	}
	x10 = x11;
	x11 = x12;
	x12 = p12 - p13;
	if (x12 < 0) {
		x12 += 2147483647;
	}

	/* Component 2 */
	h__ = x20 / 3976;
	p23 = (x20 - h__ * 3976) * 539608 - h__ * 2071;
	h__ = x22 / 24919;
	p21 = (x22 - h__ * 24919) * 86098 - h__ * 7417;
	if (p23 < 0) {
		p23 += 2145483479;
	}
	if (p21 < 0) {
		p21 += 2145483479;
	}
	x20 = x21;
	x21 = x22;
	x22 = p21 - p23;
	if (x22 < 0) {
		x22 += 2145483479;
	}

	/* Combination */
	z__ = x12 - x22;
	if (z__ <= 0) {
		z__ += 2147483647;
	}
	ret_val = z__ * 4.656612873077392578125e-10;
	return ret_val;
}
