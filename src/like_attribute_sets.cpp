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
 

#include "like_attribute_sets.hpp"

using namespace bull;

void bull::delVecLikeAttr(VecLikeAttrPtr &vla)
{
	for (VecLikeAttrPtr::iterator la = vla.begin(); la != vla.end(); ++la) {
		LikeAttr * lap = *la;
		delete lap;
	}
}
void bull::delVecTreeSimAttr(VecTreeSimAttrPtr &vla)
{
	for (VecTreeSimAttrPtr::iterator la = vla.begin(); la != vla.end(); ++la) {
		TreeSimAttr * sap = *la;
		delete sap;
	}
}

void SetOfTreeSimAttr::SetNumModsInPart(unsigned partnum, unsigned modnum)
{	
	VecTreeSimAttrPtr & vla = simAttribs[partnum];
	delVecTreeSimAttr(vla);
	vla.resize(modnum);
}

void SetOfTreeSimAttr::PrintModel(std::ostream &)
{	
	for (unsigned np=0; np < simAttribs.size(); np++) {
		VecTreeSimAttrPtr & row = simAttribs[np];
		for (unsigned nm=0; nm < row.size(); nm++) {
			row[nm]->GetModel()->PrintPAUPLsetCommand();
		}
	}
}

void SetOfTreeSimAttr::InitializeLikeAttrStatics(unsigned partnum,unsigned modn)
{
	LikeAttr::currNChar = GetNCharInEachDataSet(partnum);
	LikeAttr::currNStates = GetNStates(partnum,modn);
	LikeAttr::currModelDirty = true; //TEMPORARY
}

