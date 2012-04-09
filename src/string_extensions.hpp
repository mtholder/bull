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
 

#ifndef STRING_EXTENSIONS_H
#define STRING_EXTENSIONS_H

#include <string>
#include <vector>
#include "ncl/nxsstring.h"

std::string	   &AppendChar(std::string& s, const char c);
std::string	   &AppendDouble(std::string& s, const double d);
std::string	   &AppendInt(std::string& s, const int i);
std::string	   &AppendLong(std::string& s, const long l);

std::string	   &BlanksToUnderscores( std::string& s );
int				BreakString(std::vector<std::string> &vec, std::string phrase);
char			GetNextGraphicalNonCommentedChar(std::string fullstr, std::string::iterator &c);
char			GetNextNonCommentedChar(std::string fullstr, std::string::iterator &c);
std::string		GetTreeToken(std::string fullstr, std::string::iterator &c);
bool			IsDouble(const std::string);
bool			IsInteger(const std::string &);
void			ToUpper(std::string &);
std::string	   &UnderscoresToBlanks(std::string &);


inline std::string &AppendChar(std::string& s, const char c) {
	char cs[2];
	cs[0] = c;
	cs[1] = '\0';
	s.append(cs);
	return s;
}

inline std::string &AppendInt(std::string& s, const int i) {
	char tmp[81];
	sprintf( tmp, "%d", i );
	s.append(tmp);
	return s;
}

inline std::string &AppendLong(std::string& s, const long i) {
	char tmp[81];
	sprintf( tmp, "%ld", i );
	s.append(tmp);
	return s;
}


inline bool IsInteger(std::string& s)
{
	long n;
	return NxsString::to_long(s.c_str(), &n);
}

inline bool IsDouble(std::string s)
{	
	double n;
	return NxsString::to_double(s.c_str(), &n);
}

inline void ToUpper(std::string& s)
{	
	NxsString::to_upper(s);
}


#endif




