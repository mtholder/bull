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

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "complex.hpp"
//From John's Mr Bayes
#include <cassert>
#define PIOVER2 1.57079632679489662
complex Complex(double a, double b)
{
	complex c;
	c.re = a;
	c.im = b;
	return c;
}
void free_pscmatrix(complex **m)
// free a complex matrix allocated by psdmatrix()
// this works - I don't know why its a char *, maybe it should be a void *?

{
	//free((char *) (m[0]));
	//free((char *) (m));
	if (m)
		{if(m[0])	delete []m[0];
		delete []m;
		}
}


/*----------------------------------------------------------------------------
|
|		ComplexInvertMatrix
|
|		Invert matrix 'a' using LU-decomposition technique, 
|	storing inverse in 'a_inv'.	 Matrix 'a'
|		is destroyed.  Returns 1 if matrix is singular, 0 otherwise.
*/

int ComplexInvertMatrix(complex **a, int n, double *dwork, int *indx, 
										complex **a_inv, complex *col)
//		complex		**a;		matrix represented as vector of row pointers
//		int				n;				order of matrix
//		double		*dwork;		work vector of size n
//		int				*indx;		work vector of size n
//		complex		**a_inv;	inverse of input matrix a (matrix a is destroyed)
//	complex			*col;		for the second part, involving back subst

{
	int				rc, i, j;
	
	rc = ComplexLUDecompose(a, n, dwork, indx, (double *)NULL);
	//pprintf("ComplexLUDecompose() returned %i\n", rc);

	if (rc == 0) {
		for (j = 0; j < n; j++) {
			for (i = 0; i < n; i++)
				col[i] = Complex(0.0, 0.0);
			col[j] = Complex(1.0, 0.0);
			ComplexLUBackSubst(a, n, indx, col);
			for (i = 0; i < n; i++)
				a_inv[i][j] = col[i];
		}
	}
	return rc;
}


complex Cadd(complex a, complex b)
{
	complex c;
	c.re = a.re + b.re;
	c.im = a.im + b.im;
	return c;
}
complex RCmul(double a, complex b)
{
	complex c;
	c.re = a * b.re;
	c.im = a * b.im;
	return c;
}


complex Cmul(complex a, complex b)
{
	complex c;
	c.re = a.re * b.re - a.im * b.im;
	c.im = a.im * b.re + a.re * b.im;
	return c;
}
complex Cexp(complex a)
{
   complex c;
   c.re = exp(a.re);
   if (fabs(a.im) == 0) c.im = 0; 
   else	 { c.im = c.re*sin(a.im); c.re*=cos(a.im); }
   return (c);
}


complex Csub(complex a, complex b)
{
	complex c;
	c.re = a.re - b.re;
	c.im = a.im - b.im;
	return c;
}

complex Conj(complex a)
{
	complex c;
	c.re = a.re;
	c.im = -a.im;
	return c;
}

complex Cdiv(complex a, complex b)
{
	complex c;
	double r, den;
	if ( fabs(b.re) >= fabs(b.im) ) {
		r = b.im / b.re;
		den = b.re + r * b.im;
		c.re = (a.re + r * a.im) / den;
		c.im = (a.im - r * a.re) / den;
	} else {
		r = b.re / b.im;
		den = b.im + r * b.re;
		c.re = (a.re * r + a.im) / den;
		c.im = (a.im * r - a.re) / den;
	}
	return c;
}

double Cabs(complex a)
{
	double x, y, ans, temp;
	x = fabs(a.re);
	y = fabs(a.im);
	if (x == 0.0) 
		ans = y;
	else if (y == 0.0)
		ans = x;
	else if (x > y) {
		temp = y / x;
		ans = x * sqrt(1.0 + temp * temp);
	} else {
		temp = x / y;
		ans = y * sqrt(1.0 + temp * temp);
	}
	return ans;
}

complex Csqrt(complex a)
{
	complex c;
	double x, y, w, r;
	if ( (a.re == 0.0) && (a.im == 0.0) ) {
		c.re = 0.0;
		c.im = 0.0;
		return c;
	} else {
		x = fabs(a.re);
		y = fabs(a.im);
		if (x >= y) {
			r = y / x;
			w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r * r)));
		} else {
			r = x / y;
			w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r * r)));
		}
		if (a.re >= 0.0) {
			c.re = w;
			c.im = a.im / (2.0 * w);
		} else {
			c.im = (a.im >= 0.0) ?	w : -w;
			c.re = a.im / (2.0 * c.im);
		}
		return c;
	}
}
			


complex Clog(complex a)
{
	complex c;
	c.re = log(Cabs(a));
	if (a.re == 0.0) {
		c.im = PIOVER2;
	} else {
		c.im = atan2(a.im, a.re);
	}
	return c;
}

complex **pscmatrix(int dim)
// allocate a complex square matrix with subscript 
// range m[0..dim-1][0..dim-1]	Malloc's replaced with new by MTH
{
	int i;
	complex **m;
	/* allocate pointers to rows */
				//m=(complex **) malloc((size_t)((dim)*sizeof(complex*)));
	m=new complex *[dim];
	if (!m) {
		//printf("allocation error in pscmatrix 1.\n");
		exit(1);
	}
	// allocate rows and set pointers to them 
				//m[0]=(complex *) malloc((size_t)((dim*dim)*sizeof(complex)));
	m[0]=new complex[dim*dim];
	if (!m[0]) {
	   // printf("allocation error in pscmatrix 2.\n");
		exit(1);
	}
	for (i=1; i < dim; i++) {
		m[i]=m[i-1]+dim;
	}
	// return pointer to array of pointers to rows
	return m;
}

void dump_pscmatrix(complex **m, int dim)
{
	int		row, col;
	
	printf("{");
	for (row = 0; row < (dim - 1); row++) {
		printf("{");
		for (col = 0; col < (dim - 1); col++) {
			printf("%f + %f I, ", m[row][col].re, m[row][col].im);
			if (col == 1) printf("\n	");
		}
		printf("%f + %f I},\n", 
			m[row][dim - 1].re, m[row][dim - 1].im);
	}
	printf("{");
	for (col = 0; col < (dim - 1); col++) {
		printf("%f + %f I, ", m[dim - 1][col].re, m[dim - 1][col].im);
		if (col == 1) printf("\n	");
	}
	printf("%f + %f I}}", 
		m[dim - 1][dim - 1].re, m[dim - 1][dim - 1].im);
	printf("\n");
}

void copy_pscmatrix(complex **from, complex **to, int dim)
{
	int row, col;

	for (row = 0; row < dim; row++) {
		for (col = 0; col < dim; col++) {
			to[row][col] = from[row][col];
		}
	}

}

void copy_pscmatrix(complex *from, complex *to, int dim)
{	//MTH version with pointer math
	assert(dim < 180);
	int dimt=dim*dim;
	for (int i=0; i < dimt; i++)
		*to++=*from++;

}

void dump_complexVector(complex *vec, int dim)
{
	int		i;
	
	printf("{");
	for (i = 0; i < (dim - 1); i++) {
		printf("%f + %f I, ", vec[i].re, vec[i].im);
		if (i == 1) printf("\n	  ");
	}
	printf("%f + %f I}\n", vec[dim - 1].re, vec[dim - 1].im);
}




/*-------------------------------------------------------------------------
|
|		ComplexLUDecompose
|
|		Replace matrix 'a' with its LU-decomposition.  
|	Returns 1 if matrix is singular, 0
|		otherwise.
*/

int ComplexLUDecompose(complex **a, int n, double *vv, int *indx, double
*pd)
//		complex		**a;	the matrix whose LU-decomposition is wanted
//		int				n;			order of a
//		double		*vv;	work vector of size n (stores implicit 
//								scaling of each row)
//		int				*indx;	=> row permutation according to partial 
//								pivoting sequence
//		double		*pd;	=> 1 if number of row interchanges was even, 
//								-1 if odd (NULL OK)
{
	int				i, imax, j, k;
	double		big, dum, temp, d;
	complex		sum, cdum;

	d = 1.0;
	imax = 0; // only to shut the compiler up.

	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++) {
			if ((temp = Cabs(a[i][j])) > big)
				big = temp;
		}
		if (big == 0.0) {
				printf("singular matrix in routine ComplexLUDecompose");
			return 1;
		}
		vv[i] = 1.0 / big;
	}
	
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++) 
				sum = Csub(sum, Cmul(a[i][k], a[k][j]));
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			sum = a[i][j];
			for (k = 0; k < j; k++)
				sum = Csub(sum, Cmul(a[i][k], a[k][j]));
			a[i][j] = sum;
			dum = vv[i] * Cabs(sum);
			if (dum >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				cdum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = cdum;
			}		
			d = -d;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j].re == 0.0 && a[j][j].im == 0.0)
			a[j][j] = Complex(1.0e-20, 1.0e-20);
		if (j != n - 1){
			cdum = Cdiv(Complex(1.0, 0.0), a[j][j]);
			for (i = j + 1; i < n; i++)
				a[i][j] = Cmul(a[i][j], cdum);
		}
	}

	if (pd != NULL)
		*pd = d;
	return 0;
}

/*----------------------------------------------------------------------
|
|		ComplexLUBackSubst
|
|		Perform back-substition into LU-decomposed matrix in order 
|		to obtain inverse.
*/

void ComplexLUBackSubst(complex **a, int n, int *indx, complex *b)
{
	int				i, ip, j,
				ii = -1;
	complex		sum;

	for (i = 0; i < n; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii >= 0) {
			for (j = ii; j <= i - 1; j++)
				sum = Csub(sum, Cmul(a[i][j], b[j]));
		} else if ((sum.re != 0.0) || (sum.im != 0.0))
			ii = i;
		b[i] = sum;
	}
	
	for (i = n - 1; i >= 0; i--) {
		sum = b[i];
		for (j = i + 1; j < n; j++)
			sum = Csub(sum, Cmul(a[i][j], b[j]));
		b[i] = Cdiv(sum, a[i][i]);
	}
}


