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
 

#ifndef GENETICCODESDOTH
#define GENETICCODESDOTH


#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 
namespace bull {

class GeneticCode 
{
	int encodingType;
	char **codonOrdToAANex; //null terminated string of the amino acid abbreviation of each codon
	int *codNum; //number of the amino acid coded by this codon
	int *nCodByAA; //number of codons that code for this amino acid
	int **codByAA; //array of the codon indeces that code for this amino acid (65 means you've gone off the end the array (the safer way is to use NCodByAA array to find out how many codons there are)
	public :
	GeneticCode(const char *i,int et); //a string of single letter aa abbreviations in the left to right order of the standard table form of codons
	int *GetCodedAANum()	{return codNum; }
	int **GetCodonsByAA()	{return codByAA; }
	char **GetCodonOrdToAANex() {return codonOrdToAANex; }
	int *GetNCodByAA()	{return nCodByAA; }
	~GeneticCode();
	void FillProtParsStepMatrix(int **mat);
};

int GetCodonIndexFromStandardTableIndex(int sti);
int GetBaseIndexFromStandardBaseIndex(int sti);
int GetAminoAcidNumber(char c);

} //namespace bull

#endif

