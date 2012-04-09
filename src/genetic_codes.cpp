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
 

#include <cmath>
#include "genetic_codes.hpp"
#include "mth_exception.hpp"


using namespace bull;



GeneticCode::GeneticCode(const char *inputStr,int et)
{	encodingType=et;
	codonOrdToAANex=new char*[64];
	codonOrdToAANex[0]=new char[128];
	codonOrdToAANex[0][1]='\0';
	for (int i=1; i < 64; i++)
		{codonOrdToAANex[i]=codonOrdToAANex[i-1]+2;
		codonOrdToAANex[i][1]='\0';
		}
	codByAA=new int*[21];
	codByAA[0]=new int[126];
	for (int i=1; i < 21; i++)
		codByAA[i]=codByAA[i-1]+6;
	codNum=new int[64];
	nCodByAA=new int [21];
	
	//the input string is in the odd order of U,C,A,G  sorted first base, third base,second base
	for (int i=0; i < 64; i++)
		{int mycodon=GetCodonIndexFromStandardTableIndex(i);
		codonOrdToAANex[mycodon][0]=inputStr[i];
		codNum[mycodon]=GetAminoAcidNumber(inputStr[i]);
		}
	for (int i=0; i < 21; i++)
		{nCodByAA[i]=0;
		for (int j=0; j < 64; j++)
			{if(codNum[j] == i)
				{codByAA[i][nCodByAA[i]]=j;
				nCodByAA[i]++;
				}
			}
		if (nCodByAA[i]>6)
			throw MTHException("Too many degenerate codons in GeneticCode");
		for (int k=nCodByAA[i]; k < 6; k++)
			codByAA[i][k]=65;
		}
				
		
}
GeneticCode::~GeneticCode()
{	delete [] codonOrdToAANex[0];
	delete [] codonOrdToAANex;
	delete [] codByAA[0];
	delete [] codByAA;
	delete [] codNum;
	delete [] nCodByAA;		
}

int bull::GetCodonIndexFromStandardTableIndex(int sti)
{	int fbase=(int) floor(((double)sti+.001)/16);
	sti%=16;
	int tbase=(int)floor(((double)sti+.001)/4);
	int sbase=sti%4;
	int transFbase=GetBaseIndexFromStandardBaseIndex(fbase);
	int transSbase=GetBaseIndexFromStandardBaseIndex(sbase);
	int transTbase=GetBaseIndexFromStandardBaseIndex(tbase);
	return (16*transFbase+transSbase*4+transTbase);
}

int bull::GetBaseIndexFromStandardBaseIndex(int sti)// the order in the "standard" U C A G in codon tables
{	if (sti == 0)
		sti=3;
	else if (sti == 2)
		sti=0;
	else if (sti == 3)
		sti=2;
	return sti;
}

int bull::GetAminoAcidNumber(char c)
{	if (c >= 'a' &&c <= 'Z')
	c+='A'-'a';
	switch (c)	
		{case 'A' : return 0;
		case 'C' : return 1;
		case 'D' : return 2;
		case 'E' : return 3;
		case 'F' : return 4;
		case 'G' : return 5;
		case 'H' : return 6;
		case 'I' : return 7;
		case 'K' : return 8;
		case 'L' : return 9;
		case 'M' : return 10;
		case 'N' : return 11;
		case 'P' : return 12;
		case 'Q' : return 13;
		case 'R' : return 14;
		case 'S' : return 15;
		case 'T' : return 16;
		case 'V' : return 17;
		case 'W' : return 18;
		case 'Y' : return 19;
		case '*' : return 20;
		default : throw MTHException("Unknown amino acid single letter abbrev.");
		}
}

void GeneticCode::FillProtParsStepMatrix(int **mat) {
	//fills the 21 by 21 matrix with the min # mutations required to go from one
	//amino acid to another
	for (int i=0; i < 21; i++) {
		for (int j=0; j < 21; j++) {
			mat[i][j]=0;
		}
	}

	//fill in all one step amino acids
	for (int i=0; i < 21; i++) {
		for (int j=0; j < nCodByAA[i]; j++) {
			int scod=codByAA[i][j];
			assert(scod!=65);
			int sofbase=(int)floor(((double) scod+.001)/16.0);
			scod%=16;
			int sosecbase=(int)floor(((double) scod+.001)/4.0);
			int sothibase=scod%4;
			for (int dfbase=0; dfbase < 4; dfbase++) {
				if (dfbase!=sofbase) {
					int codonIndex=(16*dfbase+sosecbase*4+sothibase);
					if (codNum[codonIndex]!=i)
						mat[i][codNum[codonIndex]]=1;
				}
			}
	
			for (int dsbase=0; dsbase < 4; dsbase++) {
				if (dsbase!=sosecbase) {
					int codonIndex=(16*sofbase+dsbase*4+sothibase);
					if (codNum[codonIndex]!=i)
						mat[i][codNum[codonIndex]]=1;
				}
			}
	
			for (int dtbase=0; dtbase < 4; dtbase++) {
				if (dtbase!=sothibase) {
					int codonIndex=(16*sofbase+sosecbase*4+dtbase);
					if (codNum[codonIndex]!=i)
						mat[i][codNum[codonIndex]]=1;
				}
			}
		}
	}
	
	//fill in all two mutation steps
	for (int i=0; i < 21; i++) {
		for (int j=i+1; j < 21; j++) {
			if (mat[i][j] == 0) {
				bool found=false;
				for (int k=0; !found && k < 21; k++) {
					if (mat[i][k] == 1 && mat[k][j] == 1)
						found=true;
				}
				if (found) {
					mat[i][j]=2;
					mat[j][i]=2;
				}
			}
		}
	}

	for (int i=0; i < 21; i++) {
		for (int j=i+1; j < 21; j++) {
			if (mat[i][j] == 0) {
				mat[i][j] = 3;
				mat[j][i] = 3;
			}
		}
	}
 }

