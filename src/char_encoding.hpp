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

#ifndef CHARENCODING
#define CHARENCODING


#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 

#include "mth_exception.hpp"
namespace bull {

class UnknownEncoding : public MTHException 
{
	public: 
		UnknownEncoding(const char *p)
			:MTHException(p) 
			{}
};

enum EncodingType {
	DNANoGap=0 , 
	MitoCodons=1 , 
	AminoAcid=2 , 
	NucCodons
};

std::string DecodeDNACharAsStr(short c,bool keepGap/*=false*/);
int			DecodeDNACharInOneChar(short c,bool keepGap/*=false*/);
std::string DecodeProteinCharAsStr(short *c,bool keepGap/*=false*/);
int			DecodeProteinCharInOneChar(short *c,bool keepGap/*=false*/);
std::string DecodeState(short *state,int typeofEncoding);
void		FillProtParsStepMatrix(int **mat,int m);
char	  **GetOrdinationToNexusTranslator(int fromtype,int toType);
int		  **GetCodByAA(int m);
int		   *GetNCodByAA(int m);
int		   *GetCodNum(int m);
bool		IsGap(short *,int typeOfEncoding);
std::string NexusNameOfEncoding(int outputEncodingType);
int			NumDNACharactersPerCharacter(int i);
int			NumColumnsPerCharacter(int i);
int			NumShortsPerCharacter(int i);
int			NumStatesInLastShort(int i);
} // namespace bull
#endif
