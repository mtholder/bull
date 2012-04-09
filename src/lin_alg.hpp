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
 

#ifndef JHLINALG
#define JHLINALG

#include "basic_bull.hpp"

/*	lin_alg.hpp
|
|	Prototypes for matrix-inversion and eigensystem functions
|
|	Copyright (c) 1998 by David L. Swofford, Smithsonian Institution.
|	All rights reserved.
|
|	NOTE: if ANSI function prototypes are not supported, define NO_PROTOTYPES
|		  before including this file.
*/
#include "complex.hpp"
#include "matrices.hpp"

#define RC_COMPLEX_EVAL 2	/* code that complex eigenvalue obtained */

extern int	InvertMatrix (double **a, int n, double *col, int *indx, double **a_inv);
extern int	LUDecompose (double **a, int n, double *vv, int *indx, double *pd);
extern int	EigenRealGeneral (int n, double **a, double *v, double *vi, double **u, int *iwork, double *work);
int GetEigens (int n, double **qMatrix, double *eigenValues, double *eigvalsImag, double **eigvecs, double **inverseEigvecs, complex **Ceigvecs, complex **CinverseEigvecs);
//int ChangeMatrix (double t, double r, double **p, int n, double *eigenValues, double **eigvecs, double **inverseEigvecs);
//int ChangeMatrix (double t, double r, double **p, int n, double *eigenValues, double **eigvecs, double **inverseEigvecs);
//int ChangeMatrix (double rt, double **p, int n, double *eigenValues, double **eigvecs, double **inverseEigvecs);
void ChangeMatrixWithOutSharedMatrix(double rt, double *p, int n, double *eigenValues, double **eigvecs, double **inverseEigvecs,double *gexp);
void ChangeColumnWithOutSharedMatrix(double rt, double *p, int n, double *eigenValues, double **eigvecs, double **inverseEigvecs,double *gexp,int onlycol);
int ComplexChangeMatrix (int n, double **p, double rt, double *eigenValues, double *eigvalsImag, complex **Ceigvecs, complex **CinverseEigvecs);
#ifdef CONDENSE_MATRICES
	void ChangeMatrix (double rt, double *p, int n, double *eigenValues,  double *EigInvMult, double *gexp);
	void ChangeMatrix (double rt, double *p, int n, double *eigenValues,  double *EigInvMult, double *gexp,int nonzero,double *preCalced);
	void ChangeColumn (double rt, double *p, int n, double *eigenValues,  double *EigInvMult, double *gexp,int onlycol);
	void ChangeColumn (double rt, double *p, int n, double *eigenValues,  double *EigInvMult, double *gexp,int nonzero,double *preCalced,int onlycol);
	void ChangeRow (double rt, double *p, int n, double *eigenValues,  double *EigInvMult, double *gexp,int onlyrow);
	void ChangeRow (double rt, double *p, int n, double *eigenValues,  double *EigInvMult, double *gexp,int nonzero,double *preCalced,int onlyrow);
	void CalculateAndCondenseEigInvEigMult(double *eigvecs, double **inverseEigvecs,int n,double *EigInvMult);
#else
	void ChangeMatrix (double rt, double *p, int n, double *eigenValues,  double *EigInvMult, double *gexp);
	void ChangeMatrix (double rt, double *p, int n, double *eigenValues,  double *EigInvMult, double *gexp,int nonzero,double *preCalced);
	void ChangeColumn (double rt, double *p, int n, double *eigenValues,  double *EigInvMult, double *gexp,int onlycol);
	void ChangeColumn (double rt, double *p, int n, double *eigenValues,  double *EigInvMult, double *gexp,int nonzero,double *preCalced,int onlycol);
	void ChangeRow (double rt, double *p, int n, double *eigenValues,  double *EigInvMult, double *gexp,int onlyrow);
	void ChangeRow (double rt, double *p, int n, double *eigenValues,  double *EigInvMult, double *gexp,int nonzero,double *preCalced,int onlyrow);
	void CalculateGlobalEigInvEigMult(double *eigvecs, double **inverseEigvecs,int n,double *EigInvMult);
	void CalculateReorderedEigInvEigMult(double *eigvecs, double **inverseEigvecs,int n,double *EigInvMult);
	void DoPreSummation(double *EigInvMult,double *PreSum,int nlast,int n);
#endif


#endif

