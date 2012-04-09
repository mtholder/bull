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
  

#include "basic_bull.hpp"
#include "char_encoding.hpp"
#include "genetic_codes.hpp"
#include "string_extensions.hpp"



using namespace bull;


char *dnaOrdToDNANex[]={"A","C","G","T","-"};
char *aaOrdToAANex[]={"A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*"};
char *codonOrdToDNANex[]={	"AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT",
								"CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT",
								"GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT",
								"TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"};

GeneticCode MitoCode("FSYCFSYCLS*WLS*WLPHRLPHRLPQRLPQRITNSITNSMTK*MTK*VADGVADGVAEGVAEG",EncodingType(MitoCodons));
GeneticCode NucCode("FSYCFSYCLS**LS*WLPHRLPHRLPQRLPQRITNSITNSITKRMTKRVADGVADGVAEGVAEG",EncodingType(NucCodons));
int bull::DecodeDNACharInOneChar(short c,bool keepGap/*=false*/)
{	switch(c)	{  case 1 : return 'A';
			   case 2 : return 'C';
			   case 4 : return 'G';
			   case 8 : return 'T';
			   case 15 : return 'N';
			   case 16 :  if (keepGap)	return '-';
						  else			throw MTHException("Incorrect DNA Character");
			   case 5 :	  return 'R';
			   case 10 :   return 'Y';
			   case 3 :	  return 'M';
			   case 6 :	  return 'S';
			   case 7 :	  return 'V';
			   case 9 :	  return 'W';
			   case 11 :   return 'H';
			   case 12 :   return 'K';
			   case 13 :   return 'D';
			   case 14 :   return 'B';
			   case 31 :   if (keepGap) return '-';
							//else falls through to throw exception
			   default : throw MTHException("Incorrect DNA Character");
			}
	return(0);
}
int *bull::GetCodNum(int m)
{	if (m == EncodingType(MitoCodons))
		return MitoCode.GetCodedAANum();
	if (m == EncodingType(NucCodons))
		return NucCode.GetCodedAANum();
	throw MTHException("Unknown genetic code");
}


int **bull::GetCodByAA(int m)
{	if (m == EncodingType(MitoCodons))
		return MitoCode.GetCodonsByAA();
	if (m == EncodingType(NucCodons))
		return NucCode.GetCodonsByAA();
	throw MTHException("Unknown genetic code");
}

int *bull::GetNCodByAA(int m)
{	if (m == EncodingType(MitoCodons))
		return MitoCode.GetNCodByAA();
	if (m == EncodingType(NucCodons))
		return NucCode.GetNCodByAA();
	throw MTHException("Unknown genetic code");

}


char **bull::GetOrdinationToNexusTranslator(int fromType,int toType)
{	if (toType == EncodingType(DNANoGap))
		{if(fromType == EncodingType(DNANoGap))
			return dnaOrdToDNANex;
		if (fromType == EncodingType(MitoCodons) || fromType == EncodingType(NucCodons))
			return codonOrdToDNANex;
		throw MTHException("The requested OrdToNexTranslator hasn't been written yet");
		}
	if (toType == EncodingType(AminoAcid))
		{if(fromType == EncodingType(AminoAcid))
			return aaOrdToAANex;
		else if (fromType == EncodingType(MitoCodons))
				return MitoCode.GetCodonOrdToAANex();
		else if (fromType == EncodingType(NucCodons))
				return NucCode.GetCodonOrdToAANex();
		throw MTHException("The requested OrdToNexTranslator hasn't been written yet");
		}
	throw MTHException("The requested OrdToNexTranslator hasn't been written yet");
}
std::string bull::DecodeDNACharAsStr(short c,bool keepGap/*=false*/)
{	std::string retstr;
	retstr=(char)  DecodeDNACharInOneChar(c,keepGap);
	return retstr;
}

std::string bull::DecodeProteinCharAsStr(short *c, bool keepGap/*=false*/) {
	try {
		const int i = DecodeProteinCharInOneChar(c, keepGap);
		std::string s;
		AppendInt(s, i);
		return s;
	}
	catch (MTHException){}//not just one character or a standard one character ambiguity code
	short smask=1,temp[2];
	temp[1]=0;
	std::string retstr="(";
	for (int i=0; i < 16; i++)
		{if(*c&smask)
			{temp[0]=smask;
			retstr+=(char) DecodeProteinCharInOneChar(temp,keepGap);
			}
		smask <<= 1;
		}
	smask=1;
	c++;
	temp[0]=0;
	for (int i=0; i < 6; i++)
		{if(*c&smask)
			{temp[1]=smask;
			retstr+=(char) DecodeProteinCharInOneChar(temp,keepGap);
			}
		smask <<= 1;
		}
	retstr+=")";
	return retstr;
}

int bull::DecodeProteinCharInOneChar(short *c,bool keepGap/*=false*/)
{	short maxs=16384;
	maxs <<= 1; //maxs ==  1000 0000 0000 0000
	if (c[1] == 0)
		switch(*c)	{case 1 :	return 'A';
					case 2 :	return 'C';
					case 4 :	return 'D';
					case 8 :	return 'E' ;
					case 16 :	return 'F' ;
					case 32 :	return 'G' ;
					case 64 :	return 'H' ;
					case 128 :	return 'I' ;
					case 256 :	return 'K' ;
					case 512 :	return 'L' ;
					case 1024 : return 'M' ;
					case 2048 : return 'N' ;
					case 4096 : return 'P';
					case 8192 : return 'Q';
					case 16384 :	return 'R' ;
					default :
								if (*c == maxs)	
									return 'S';
								throw MTHException("Incorrect Protein Character");
					}
	if ((*c) == 0)
		switch(c[1])  
				{case 1 :	return 'T';
				case 2 :	return 'V';
				case 4 :	return 'W'; 
				case 8 :	return 'Y'; 
				case 16 :	return '*'; 
				case 32 : if (keepGap)	return '-';
						 //else is the fall through exception
				default : throw MTHException("Incorrect Protein Character");
				}
			  
   if (c[0] == ~0)
			{if(keepGap && c[1] == 31)	return 'X';
			if (!keepGap && c[1] == 63)	return '?';
			throw MTHException("Incorrect Protein Character");
			}
   if (c[0] == 2052 && c[1] == 0) return 'B';
   if (c[0] == 8200 && c[1] == 0) return 'Z';
	throw MTHException("Incorrect Protein Character");
}


std::string bull::NexusNameOfEncoding(int i)
{	std::string tn;
	if (i == EncodingType(DNANoGap) || i == EncodingType(MitoCodons) || i == EncodingType(NucCodons)) 
			tn="DNA";
	if (i == EncodingType(AminoAcid)) 
		tn="PROTEIN";
	return tn;
}

int bull::NumDNACharactersPerCharacter(int i)	{
	if (i == EncodingType(DNANoGap))	return 1;
	if (i == EncodingType(MitoCodons) || i == EncodingType(NucCodons)|| i == EncodingType(AminoAcid)) return 3;
	throw UnknownEncoding("Not DNANogap Or Codons Or Amino Acid encoding");
	return 0;
}

int bull::NumColumnsPerCharacter(int i) {
	if (i == EncodingType(DNANoGap) || i == EncodingType(AminoAcid)) return 1;
	if (i == EncodingType(MitoCodons) || i == EncodingType(NucCodons))	return 3;
	throw UnknownEncoding("Not DNANogap Or Codons Or Amino Acid encoding");
	return 0;
}

std::string bull::DecodeState(short *inp,int datatype)
{	if (datatype == EncodingType(AminoAcid))		return DecodeProteinCharAsStr(inp,true);
	else if (datatype == EncodingType(DNANoGap))	return DecodeDNACharAsStr(*inp,true);
	throw	MTHException("Entered unwritten code DecodeState" );
}

