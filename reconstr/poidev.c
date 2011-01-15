/* The poidev routine gives a poisson-distributed random value.
   aruments : float xm, long * idum
   xm : mean expectation
   idum : pointer to negative integer (always the same pointer)
*/

/*   This file is part of a software package written by 
     Rainer Heintzmann
     Institute of Applied Optics and Information Processing
     Albert Ueberle Strasse 3-5
     69120 Heidelberg
     Tel.: ++49 (0) 6221  549264
     No garantee, whatsoever is given about functionallaty and  savety.
     No warranty is taken for any kind of damage it may cause.
     No support for it is provided !

     THIS IS NOT FREE SOFTWARE. All rights are reserved, and the software is only
     given to a limited number of persons for evaluation purpose !
     Please do not modify and/or redistribute this file or any other files, libraries
     object files or executables of this software package !


    The code in this module was directly obtained from the Numerical Recipies in C book.
    The copyright stated there does therefore also apply to this code.
*/


#include "poidev.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846L 
#endif

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/** generates random numbers from 0.0 to 1.0 */
float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = - (*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum = IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j]= * idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

float gammln(float xx)
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

float poidev(float xm, long * idum)
{
static float sq, alxm,g,oldm=(-1.0);
float em,t,y;

if (xm < 12.0) {
	if (xm != oldm) {
		oldm=xm;
		g=exp(-xm);
	}
	em = -1;
	t=1.0;
	do {
		++em;
		t *= ran1(idum);
	} while (t > g);
} else {
	if (xm != oldm) {
		oldm=xm;
		sq=sqrt(2.0*xm);
		alxm=log(xm);
		g=xm*alxm-gammln(xm+1.0);
	}
	do {
		do {
			y=tan(M_PI*ran1(idum));
			em=sq*y+xm;
		} while (em < 0.0);
		em=floor(em);
		t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
	} while (ran1(idum) > t);
}
return em;
}

/*
#define MEAN   3.0
#define MAX   100
#define RANDS 100000


int main(void)
{
int i,ival;
long rval=-5;

int array[MAX];
float val=0.0;

for (i=0;i<MAX;i++)
	array[i]=0;

for (i=0;i<RANDS;i++)
	{
	val = poidev(MEAN, & rval);
	ival = (int) val;
	// printf(" %g %d\n",val,ival);
	if (ival < MAX) array[ival] = array[ival] + 1;
	}

for (i=0;i<MAX;i++)
    printf(" %g \n",((double) array[i])/(double) RANDS);

return;
}

*/

