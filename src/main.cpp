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
 

#include<iostream>
#include<iomanip>
#include <ctime>
#include "bull.hpp"
#include "tools.hpp"
#include "ssrf_codon_sub_mod.hpp"

using namespace std;

void SetDefaults(void);
long bull::gRngSeed;

int main(int argc, char* argv[]) {

	cout << setprecision(8);
	char *infile = NULL;
	if (argc > 2) {
		cerr << "Sorry, this program can accept at most one command" << endl;
		cerr << "line argument, which must be the name of a NEXUS" << endl;
		cerr << "data file." << endl;
		cerr << endl;
		return 0;
		}
	else if (argc > 1) {
		infile = argv[1];
	}
	SetDefaults();
	
	bull::BullKernel kernel;
	bull::BullShell bull(kernel);
	try {
		bull.Run(infile);
	}
	catch(bull::MTHException x) {
		cerr<< x.msg << endl;
		return 1;
	}
	catch(...) {
		cerr << "Exception of undetermined type caught." << endl;
		return 1;
	}						
	return 0;
}


void SetDefaults(void){
	bull::RateManager::set_default(0.1,0.5);
	bull::gRngSeed=(long) time(NULL);
	for (int i=0; i < 10; i++)
		RandomNumber(); 
}


