// This file is part of BULL, a program for phylogenetic simulations
// most of the code was written by Mark T. Holder.

//	This program is for internal use by the lab of Dr. Tandy Warnow only.
//	Do not redistribute the code.  It is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
//
//	Some of the code is from publically available source by Paul Lewis, Ziheng Yang, 
//	John Huelsenbeck, David Swofford , and others  (as noted in the code).
//	In fact the main structure of the program was created by modifying Paul Lewis' 
//	basiccmdline.cpp from his NCL
//
//	This code was used in Mark's dissertation, some changes were made in order to 
//	get it to compile on gcc.  It is possible that this porting introduced bugs (very
//	little debugging has been done on UNIX platforms).	I would suggest checking 
//	the simulator by generating data on trees with short branches, etc.

#include <iostream>
#include "util.hpp"

using std::cout;
using std::endl;
using namespace bull;

#ifndef BASIC_FUNCS
#include "util.hpp"
#endif
int itmax = 100;
double epsilon = 3.0e-7;

double gammln( double xx )
{
	double x,tmp,ser;
	static double cof[6]={76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;

	x=xx-1.0;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.0;
	for (j=0; j < 6; j++) {
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}

//
// from the Numerical Recipes (Press et al.) function of the same name
//
void gcf( double& gammcf, double a, double x, double& gln )
{
	int n;
	double gold=0.0,g,fac=1.0,b1=1.0;
	double b0=0.0,anf,ana,an,a1,a0=1.0;

	gln=gammln(a);
	a1=x;
	for (n=1; n<=itmax; n++) {
		an=(double) n;
		ana=an-a;
		a0=(a1+a0*ana)*fac;
		b0=(b1+b0*ana)*fac;
		anf=an*fac;
		a1=x*a0+anf*a1;
		b1=x*b0+anf*b1;
		if (a1) {
			fac=1.0/a1;
			g=b1*fac;
			if (fabs((g-gold)/g) < epsilon) {
				gammcf = exp( -x + a * log(x) - gln ) * g;
				return;
			}
			gold=g;
		}
	}
	throw MTHException("Error in routine GCF: a too large, itmax too small");
}


//
// from the Numerical Recipes (Press et al.) function of the same name
//
void gser( double& gamser, double a, double x, double& gln )
{
	int n;
	double sum,del,ap;

	gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) {
			cout << "\nError in routine GSER: x less than 0" << endl;
			cout << "Program aborted." << endl;
			exit(1);
	  }
		gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1; n<=itmax; n++) {
			ap += 1.0;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*epsilon) {
				gamser = sum * exp( -x + a * log(x) - gln );
				return;
			}
		}
		cout << endl << "Error in routine GSER:	 a too large, itmax too small" << endl;
		cout << "Program aborted." << endl;
		exit(1);
	}
}

//
// from the Numerical Recipes (Press et al.) function of the same name
//
double gammq( double a, double x )
{
	double gamser, gammcf, gln;

	if (x < 0.0 || a <= 0.0) {
		cout << endl << "Error in routine GAMMQ:  Invalid arguments" << endl;
		cout << "Program aborted." << endl;
		exit(1);
	}
	if (x < (a+1.0)) {
		gser(gamser, a, x, gln);
		return 1.0 - gamser;
	} else {
		gcf(gammcf, a, x, gln);
		return gammcf;
	}
}
int GetRandomIndexFromFreqs(int maxInd, const double * const *freqPP)//the assumption is that the freqs add up to 1, the maxInd is returned if the sum is lower
{	double x = RandomNumber(); 
	int i;
	for (i=0; i < maxInd-1; i++)
		{x-=**freqPP++;
		if (x < 0.0)
			return i;
		}
	return i;
}

int GetRandomIndexFromFreqs(int maxInd, const double * freqP)//the assumption is that the freqs add up to 1, the maxInd is returned if the sum is lower
{	double x=RandomNumber();
	int i;
	for (i=0; i < maxInd-1; i++)
		{x-=*freqP++;
		if (x < 0.0)
			return i;
		}
	return i;

}	


