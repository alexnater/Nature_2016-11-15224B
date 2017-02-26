/*
 * TRandomGenerator.cpp
 *
 *  Created on: Sep 24, 2009
 *      Author: wegmannd
 */

#include "TRandomGenerator.h"

//---------------------------------------------------------------------------------------
long TRandomGenerator::get_randomSeedFromCurrentTime(long* addToSeed){
   time_t seconds;
   seconds = time (NULL);
   seconds = seconds - 38*365*24*3600; //substract 38 years...
   seconds = (float) (*addToSeed+seconds);
   return seconds ;
}


#define MBIG 1000000000L
#define MSEED 161803398L
#define MZ 0
#define FAC (1.0/MBIG)

double TRandomGenerator::ran3(long *idum)
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if ((*idum < 0) || (iff == 0) ) {
                iff=1;
                mj=MSEED-(*idum < 0 ? -*idum : *idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;++i) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k<=4;++k)
                        for (i=1;i<=55;++i) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                *idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
//-------------------------------------------------------
/* binomial distribution */
float TRandomGenerator::getBiomialRand(float pp, int n){
 int j;
 static int nold=(-1);
 float am,em,g,angle,p,bnl,sq,t,y;
 static float pold=(-1.0),pc,plog,pclog,en,oldg;

 p=(pp <= 0.5 ? pp : 1.0-pp);
 am=n*p;
 if (n < 25) {
  bnl=0;
  for (j=1;j<=n;j++)
   if (ran3(&_Idum) < p) ++bnl;
 } else if (am < 1.0) {
  g=exp(-am);
  t=1.0;
  for (j=0;j<=n;j++) {
   t *= ran3(&_Idum);
   if (t < g) break;
  }
  bnl=(j <= n ? j : n);
 } else {
  if (n != nold) {
   en=n;
   oldg=gammln(en+1.0);
   nold=n;
  } if (p != pold) {
   pc=1.0-p;
   plog=log(p);
   pclog=log(pc);
   pold=p;
  }
  sq=sqrt(2.0*am*pc);
  do {
   do {
    angle=3.141592654*ran3(&_Idum);
    y=tan(angle);
    em=sq*y+am;
   } while (em < 0.0 || em >= (en+1.0));
   em=floor(em);
   t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
    -gammln(en-em+1.0)+em*plog+(en-em)*pclog);
  } while (ran3(&_Idum) > t);
  bnl=em;
 }
 if (p != pp) bnl=n-bnl;
 return bnl;
}
//--------------------------------------------------------
float TRandomGenerator::gammln(float xx)
{
 double x,y,tmp,ser;
 static double cof[6]={76.18009172947146,-86.50532032941677,
  24.01409824083091,-1.231739572450155,
  0.1208650973866179e-2,-0.5395239384953e-5};
 int j;

 y=x=xx;
 tmp=x+5.5;
 tmp -= (x+0.5)*log(tmp);
 ser=1.000000000190015;
 for (j=0;j<=5;j++) ser += cof[j]/++y;
 return -tmp+log(2.5066282746310005*ser/x);
}
//--------------------------------------------------------
double TRandomGenerator::getRand(double min, double max){
	//return a random number between min and max
	 return (ran3(&_Idum) * (max - min) + min);
}
int TRandomGenerator::getRand(int min, int maxPlusOne){
	//return an random integer between min and maxPlusOne-1
	float r=1;
	while(r==1) r=getRand(); //we have a number in [0,1[
	return min+floor(r*(maxPlusOne-min));
}

long TRandomGenerator::getRand(long min, long maxPlusOne){
	//return an random integer between min and maxPlusOne-1
	double r=1.0;
	while(r==1.0) r=getRand(); //we have a number in [0,1[
	return min+floor(r*(maxPlusOne-min));
}

/* Returns a Normal random variate based on a unit variate,
   using a random generator as a source of uniform deviates.
   Adapted from the algorithm described in the book Numerical Recipes by
   Press et al.
*/
double TRandomGenerator::getNormalRandom (double dMean, double dStdDev){
  double w, x1, x2;
   do {
      x1 = 2. * ran3(&_Idum) - 1.;
      x2 = 2. * ran3(&_Idum) - 1.;
      w = x1*x1 + x2*x2;
   } while (w >= 1. || w < 1E-30);

   w = sqrt((-2.*log(w))/w);
   x1 *= w;
   return (x1 * dStdDev + dMean);
} /* NormalRandom */
//------------------------------------------------------------------------------
/*
Rogers comments:

In answer to your question about my algorithm, I'm going to append my
entire gamma_dev function.  The comments at the top provide references
to the original sources.  The algorithm I use is supposedly the most
commonly used when alpha<1.

In case it is relevant, let me tell you about some of the trouble I've
run into generating gamma deviates with small values of alpha.  My
first gamma_dev function was in single precision.  It behaved very
strangely.  When alpha<0.1, the number of segregating sites went *up*
as alpha went *down*, which makes no sense at all.  I couldn't find
any error in the code, but I noticed that the code does things that
may stretch the limits of floating point arithmetic.  So I recompiled
using double precision for all variables within gamma_dev.  The
strange behavior went away.

The literature doesn't say much about the stability of these
algorithms when alpha is very small.  It seems that no one has ever
been interested in that case.  I'll bet that none of the commercial
statistical packages have tested their gamma deviate generator with
very small alpha values either.  Consequently, we can't test our
algorithms by comparing the quantiles of our generated values with
those generated by, say, SPSS.  The only sure way is to calculate
quantiles by direct integration of the density function.  I have done
this for alpha=0.1 and am about to compare the quantiles of my numbers
with these values.  I'll let you know what happens.

Alan

PS  Here's the code along with references.  */

/****************************************************************
Random deviates from standard gamma distribution with density
         a-1
        x    exp[ -x ]
f(x) = ----------------
         Gamma[a]

where a is the shape parameter.  The algorithm for integer a comes
from numerical recipes, 2nd edition, pp 292-293.  The algorithm for
a<1 uses code from p 213 of Statistical Computing, by Kennedy and
Gentle, 1980 edition.  This algorithm was originally published in:

Ahrens, J.H. and U. Dieter (1974), "Computer methods for sampling from
Gamma, Beta, Poisson, and Binomial Distributions".  COMPUTING
12:223-246.

The mean and variance of these values are both supposed to equal a.
My tests indicate that they do.

This algorithm has problems when a is small.  In single precision, the
problem  arises when a<0.1, roughly.  That is why I have declared
everything as double below.  Trouble is, I still don't know how small
a can be without causing trouble.  Mean and variance are ok at least
down to a=0.01.  f(x) doesn't seem to have a series expansion around
x=0.
****************************************************************/
double TRandomGenerator::getGammaRand(double a, double b) {
   if(b <= 0) {
	   throw TException("Negative value received in getGammaRand()!", _FATAL_ERROR);
   }
   return getGammaRand(a)/b;
}


double TRandomGenerator::getGammaRand(double a) {

  int ia;
  double u, b, p, x, y=0.0, recip_a;

  if(a <= 0) {
	  throw TException("Negative value received in getGammaRand()!", _FATAL_ERROR);
  }

  ia = (int) floor(a);  /* integer part */
  a -= ia;        /* fractional part */
  if(ia > 0) {
    y = getGammaRand(ia);  /* gamma deviate w/ integer argument ia */
    if(a==0.0) return(y);
  }

  /* get gamma deviate with fractional argument "a" */
  b = (M_E + a)/M_E;
  recip_a = 1.0/a;
  for(;;) {
    u = ran3(&_Idum);
    p = b*u;
    if(p > 1) {
      x = -log((b-p)/a);
      if( ran3(&_Idum) > pow(x, (double) a-1.0)) continue;
      break;
    }
    else {
      x = pow(p, recip_a);
      if( ran3(&_Idum) > exp(-x)) continue;
      break;
    }
  }
  return(x+y);
}

double TRandomGenerator::getGammaRand(int ia){
	/****************************************************************
	gamma deviate for integer shape argument.  Code modified from pp
	292-293 of Numerical Recipes in C, 2nd edition.
	****************************************************************/
  int j;
  double am,e,s,v1,v2,x,y;

  if (ia < 1) throw TException("Argument below 1 in getGammaRand()!", _FATAL_ERROR);
  if (ia < 6){
    x=1.0;
    for (j=0; j<ia; j++)
      x *= ran3(&_Idum);
    x = -log(x);
  } else {
    do {
      do {
	    do{                         /* next 4 lines are equivalent */
	    	v1=2.0*ran3(&_Idum)-1.0;       /* to y = tan(Pi * uni()).     */
	    	v2=2.0*ran3(&_Idum)-1.0;
	    } while (v1*v1+v2*v2 > 1.0);
	    y=v2/v1;
	    am=ia-1;
	    s=sqrt(2.0*am+1.0);
	    x=s*y+am;
      } while (x <= 0.0);
      e=(1.0+y*y)*exp(am*log(x/am)-s*y);
    } while (ran3(&_Idum) > e);
  }
  return(x);
}

/* -----------------------------------------------------------------------------
   BetaRandom

   returns a variate that is Beta distributed on the interval [a,b]
   with shape parameters alpha and beta.

   The Beta function has two shaping parameters, alpha and beta.
   Setting these parameters to 1.5 and 1.5 yields a normal-like
   distribution, but without tails. If alpha and beta are equal to
   1 it is a uniform distribution.

   If alpha and beta are less than 1, use a rejection algorithm;
   Otherwise use the fact that if x is distributed Gamma(alpha) and y
   Gamma(beta) then x/(x+y) is Beta(alpha, beta).

   The rejection algorithm first a Beta variate is found over the
   interval [0, 1] with not the most efficient algorithm.  This is then
   scaled at the end to desired range.

   It may be tempting to re-use the second number drawn as the first
   random number of the next iteration, and simply draw one more.
   *** Don't do it.  You will produce an incorrect distribution.  You
   must draw two new numbers for the rejection sampling to be correct.

   References:
   - Ripley, Stochastic Simulations, John Wiley and Sons, 1987, p 90.
   - J.H.Maindonald, Statistical Computation, John Wiley and Sons,
     1984, p 370.
*/
double TRandomGenerator::getBetaRandom (double alpha, double beta, double a, double b){
  if (b <= a) throw TException("Bad shape or range for a beta variate!", _FATAL_ERROR);
  return (a + getBetaRandom(alpha, beta) * (b-a));   /* Scale to interval [a, b] */
}


double TRandomGenerator::getBetaRandom (double alpha, double beta) {
  double u1, u2, w;

  if (alpha <= 0 || beta <= 0) throw TException("Bad shape or range for a beta variate!", _FATAL_ERROR);

  if ((alpha < 1) && (beta < 1))
    /* use rejection */
    do {
      u1 = ran3(&_Idum); /* Draw two numbers */
      u2 = ran3(&_Idum);

      u1 = pow(u1, (double) 1.0/alpha); /* alpha and beta are > 0 */
      u2 = pow(u2, (double) 1.0/beta);

      w = u1 + u2;

    } while (w > 1.0);

  else {
    /* use relation to Gamma */
    u1 = getGammaRand(alpha);
    u2 = getGammaRand(beta);
    w  = u1 + u2;
  }

  return u1/w;

} /* BetaRandom */
