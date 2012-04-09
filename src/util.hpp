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
 

#ifndef _BASICFUNCS
#define _BASICFUNCS



#include <fstream>
#include <iomanip>
#include <list>
#include <vector>
#include <string>
#include <iterator>
#include <cstdio>

#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 


#include "tools.hpp"
#include "mth_exception.hpp"
#include "basic_bull.hpp"

/*Combo of my funcs and stuff taken from 
PAML and Mr Bayes
*/
namespace bull {

class  VariableNotDetermined : public MTHException 
{
};

} // namespace bull;



int ClassifyChange(short int brstart,short int brend,short int *plink); //Rogers and Swofford's Fixed and Potentially Same
int ClassifyChangeTerm(short int brstart,short int brend,short int *plink); //Rogers and Swofford's Fixed and Potentially Same
int ConvertCodeToIndex(short int c);
int ConvertPermLinkToChange(short *plink);
short int DecodeAChar(short int n);
int DecodeAmbiguityCode(char *dest,short int c);
void gcf( double& gammcf, double a, double x, double& gln );
double gammln( double xx );
double gammq( double a, double x );
void gser( double& gamser, double a, double x, double& gln );
double GetNumberFromNStr(std::string::iterator &c);
void GetNStrToPunc(std::string::iterator &c,char *dest);
std::string GetNStrToPunc(std::string::iterator &inc);
long NumberOfLongsLeftInMemory();
long NumberOfCharsLeftInMemory();
long NumberOfBoolsLeftInMemory();
int NextNexusWord(std::ifstream &StreamRef,char *temp);
char nexusget(std::ifstream &inp);
double point_normal( double pr );
double point_chi2 ( double pr, double v );
void PrintPMat(double *pmat);
int ProductOfSquareMatrices(double *prod, int dimen, double *first, double *second);
double RandomNumber(long int *seed);
int GetRandomIndexFromFreqs(int maxInd, const double * const * freqPP);
int GetRandomIndexFromFreqs(int maxInd, const double * freqP);
int ScreenStrForStr(std::ifstream &StreamRef,char *source,char *fragment);
int ScreenStrForStr(char *source,char *fragment);
int StateCount(short int c);
void itos(char *dest,int x);
void ltos(char *dest,long x);
void dtos(char *dest,double x);
int NextWord(std::ifstream &StreamRef,char *temp);
int UtoCstrcpy(char *ctemp,unsigned char *temp);
int CtoUstrcpy(unsigned char *ctemp, char *temp);
int CtoUPascalstrcpy(unsigned char *ctemp, char *temp);

bool StrForDoubleInRangeNoSci(char *str,double minim,double maxim);
bool StrForDoubleInRange(char *str,double minim,double maxim);
double GetPartNumberFromStr(char *source);
bool CouldTurnIntoADouble(char *str);
bool CouldTurnIntoADoubleNoSci(char *str);
bool HasADigit(char *str);
void SortDoubleArray(double *x,int len); 
char *AllocateCharAndCopyNStr(std::string instr,int len);
std::string GetNextGraphStr(std::string::iterator &nextlet, std::string::iterator endp);
short EncodeTripletToCodon(short *str);
short ConverNumAAToLetAA(short str);
bool IsGap(short c);
short GetMissing();
short GetNStatesIncludingMissing();
short GetNStatesIncludingMissingAndGaps();
void GaussJordanElimination(double **mat,double *vars,int dim);
int StateCount( int c);
int EncodeTripletToCodon(int *str);
bool PointsAreFarEnoughApart(double a,double b);

#endif

