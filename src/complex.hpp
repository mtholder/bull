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

#ifndef JHCOMPLEX
#define JHCOMPLEX
#include <cmath>

struct complex { double re;	 double im; };

typedef struct complex complex;

extern complex	   Complex (double a, double b);
extern complex	   Cadd (complex a, complex b);
extern complex	   Csub (complex a, complex b);
extern complex	   Cmul (complex a, complex b);
extern complex	   Conj (complex a);
extern complex	   Cdiv (complex a, complex b);
extern double	   Cabs (complex a);
extern complex	   Csqrt (complex a);
extern complex	   RCmul (double a, complex b);
extern complex	   Cexp (complex a);
extern complex	   Clog (complex a);
extern complex **pscmatrix (int dim);	/* mallocs a square complex matrix*/
extern void		   free_pscmatrix (complex **m); /* frees the square complex matrix*/
extern void		   dump_pscmatrix (complex **m, int dim); /* prints a square complex matrix*/
extern void		   copy_pscmatrix (complex **from, complex **to, int dim); /* copies*/
void				copy_pscmatrix(complex *from, complex *to, int dim);
extern void		   dump_complexVector (complex *vec, int dim);
extern int		   ComplexInvertMatrix (complex **a, int n, double *dwork, int *indx, complex **a_inv, complex *col);
extern int		   ComplexLUDecompose (complex **a, int n, double *vv, int *indx, double *pd);
extern void		   ComplexLUBackSubst (complex **a, int n, int *indx, complex *b);
#endif


