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
 

#ifndef SETSOFLIKEATTRIBUTES
#define SETSOFLIKEATTRIBUTES


#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 

#include "like_attributes.hpp"
#include "part_mod_index.hpp"
namespace bull {

typedef std::vector<LikeAttr*> VecLikeAttrPtr;
typedef std::vector<VecLikeAttrPtr> VecVecLikeAttrPtr;
typedef std::vector<TreeSimAttr*> VecTreeSimAttrPtr;
typedef std::vector<VecTreeSimAttrPtr> VecVecTreeSimAttrPtr;

void delVecLikeAttr(VecLikeAttrPtr &vla);
void delVecTreeSimAttr(VecTreeSimAttrPtr &vla);
class SetOfLikeAttr
{
	public:
		
		SetOfLikeAttr(unsigned nparts)
			:likeAttribs(nparts)
			{}
			
		~SetOfLikeAttr()
		{	
			for (VecVecLikeAttrPtr::iterator vla = likeAttribs.begin(); vla != likeAttribs.end(); ++vla) {
				delVecLikeAttr(*vla);
			}
		}


		Model * GetModel(unsigned p, unsigned m) const
		{
			return likeAttribs[p][m]->GetModel();
		}

		void SetNumModsInPart(unsigned partnum, unsigned modnum)
		{
			VecLikeAttrPtr &vla = likeAttribs[partnum];
			delVecLikeAttr(vla);
			vla.resize(modnum);
		}

		void AddLikeAtt(unsigned partnum, unsigned m, LikeAttr *la)
		{
			delete likeAttribs.at(partnum).at(m);
			likeAttribs[partnum][m] = la;
		}

		unsigned GetNParts() const
		{
			return likeAttribs.size();
		}

		unsigned GetNModels(unsigned partNum) const
		{
			return likeAttribs.at(partNum).size();
		}

		LikeAttr *GetLikeAtt(unsigned partnum,unsigned m)
		{
			return likeAttribs.at(partnum).at(m);
		}
	protected:
		VecVecLikeAttrPtr likeAttribs;

};

class SetOfTreeSimAttr
{
	public:
		SetOfTreeSimAttr(unsigned nparts)
			:nChar(nparts, 0),
			simAttribs(nparts)
			{}
		~SetOfTreeSimAttr()
		{	
			for (VecVecTreeSimAttrPtr::iterator vla = simAttribs.begin(); vla != simAttribs.end(); ++vla) {
				delVecTreeSimAttr(*vla);
			}
		}
		

		
		void SetNumModsInPart(unsigned partnum,unsigned modnum);
		void AddAtt(unsigned partnum, unsigned m,TreeSimAttr *tsa)
		{
			delete simAttribs.at(partnum).at(m);
			simAttribs[partnum][m] = tsa;
		}

		unsigned GetNParts()
		{
			return simAttribs.size();
		}

		unsigned GetNModels(unsigned i)
		{
			return simAttribs.at(i).size();
		}

		const Model *GetModel(unsigned p,unsigned m)
		{
			TreeSimAttr * tsa = GetSimAtt(p,m);
			return (tsa == NULL ? NULL : tsa->GetModel());
		}
		
		TreeSimAttr *GetSimAtt(unsigned p,unsigned m)
		{
			return simAttribs[p][m];
		}
		
		void InitializeLikeAttrStatics(unsigned partnum,unsigned modn);
		
		unsigned GetNCharInSimAttr(unsigned p,unsigned m)
		{
			TreeSimAttr * tsa = GetSimAtt(p,m);
			return (tsa == NULL ? 0 : tsa->GetNChar());
		}
		
		unsigned GetNCharInEachDataSet(unsigned p)
		{
			return nChar[p];
		}

		void SetNCharInEachDataSet(unsigned p,unsigned nc)
		{
		  	nChar[p] = nc;
		}
	
		unsigned GetNStates(unsigned p,unsigned m) {
			TreeSimAttr * tsa = GetSimAtt(p,m);
			return (tsa == NULL ? 0 : tsa->GetNStates());
		}
		void PrintModel(std::ostream &usrstream);
		bool GetOwnsBrLengths(bull::PartModIndex pmi)	{
			TreeSimAttr * tsa = GetSimAtt(pmi.partNum, pmi.modNum);
			return (tsa == NULL ? false : tsa->GetOwnsBrLen());
		}

	protected:
		std::vector<unsigned> nChar;

		VecVecTreeSimAttrPtr simAttribs;

	friend class Tree;
};
} // namespace bull
#endif
