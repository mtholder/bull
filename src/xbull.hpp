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

#ifndef XNEXUS_H
#define XNEXUS_H

#include <iostream>

#include "ncl/nxsexception.h"

class XBull : public NxsException
{
	public:
		XBull(const std::string & s, file_pos fp = 0, long fl = 0L, long fc = 0L)
			:NxsException(s, fp, fl, fc)
			{}
		XBull(const std::string &s, const NxsToken &t)
			:NxsException(s, t)
			{}
		XBull(const std::string &s, const ProcessedNxsToken &t)
			:NxsException(s, t)
			{}
};

#endif

