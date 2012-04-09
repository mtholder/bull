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
 

#ifndef MTHEXCEPTION
#define MTHEXCEPTION

#include <string>

#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 

namespace bull {

class  MTHException {
	public:
	std::string msg;
	MTHException()
	{
		msg=" ";
	}
	MTHException(const char *m)
	{
		msg=m;
	}
	MTHException(std::string n, const char *m)
	{
		msg=n;
		msg+=m;
	}
	
};
class CantMove : public MTHException	{
};
class LikePlateau : public MTHException {
};
class CloseToMaxDistException : public MTHException {
};
class AtPeakLike : public MTHException	{
};
class TaxonWithoutStateException : public MTHException{
	public: TaxonWithoutStateException(const char *m){
		msg=m;
		}
};

class MaximizationError : public MTHException	{
};
} // namespace bull
#endif
