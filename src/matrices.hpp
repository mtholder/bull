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
 

#ifndef JHMATRICES
#define JHMATRICES


#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 

/*	matrices.h
|
|	Header for routines defined in matrices.c
|
|	Copyright (c) 1998 by David L. Swofford, Smithsonian Institution.
|	All rights reserved.
|
|	Things that (may) need to be done before including this file:
|	= == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==
|	  - if ANSI function prototypes are not supported, define NO_PROTOTYPES
|	  - typedef type VoidPtr as void * if possible, otherwise as char *
*/

typedef void * VoidPtr;

int			AllocMatrix (VoidPtr pA, size_t elSize, unsigned rows, unsigned cols);
void			DeallocMatrix (VoidPtr pA);
short **		ShortMatrix (unsigned rows, unsigned cols, short **a, short *buffer);
double **	DoubleMatrix (unsigned rows, unsigned cols, double **a, double *buffer);
void			SetIdentityMatrix (double **a, unsigned n);
void			CopyDoubleMatrix (double **a, double **b, unsigned m, unsigned n);
void			ListDoubleMatrix (char *title, double **a, unsigned m, unsigned n);
void			ListDoubleVector (char *title, double *a, unsigned n);

double		**psdmatrix (unsigned dim);
double		***psdmatrices (unsigned dim,unsigned num);
void			free_psdmatrix (double **m);
void			free_psdmatrices (double ***m,unsigned num);
void			copy_psdmatrix (double **from, double **to, unsigned dim);
void			copy_psdmatrix (double *from, double *to, unsigned dim);

double **new_RectMats(unsigned rows,unsigned cols);
void free_RectMats(double **arr);

#endif
