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
 
#ifndef _TREE_
#define _TREE_

#include <map>
#include <algorithm>
#include <ostream>

class NxsTaxaBlockAPI;

#include "node.hpp"
#include "basic_bull.hpp"

#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 

#include "ssrf_codon_sub_mod.hpp"
#include "like_attribute_sets.hpp"

namespace bull {

class IllegalBranchLength : public MTHException { public :
	IllegalBranchLength() : MTHException() {}
	IllegalBranchLength(const char *c) :MTHException(c) { }
	};
	
class TreeEncodingException : public MTHException { public :
	TreeEncodingException() : MTHException() {}
	TreeEncodingException(const char *c) :MTHException(c) { }
	};

class Tree
{
	public:
		static unsigned timesscored;
		unsigned tempOneBrentScored;
		double likelihood;
		bool hasBranchLengths;
		Tree();
		Tree(Node *);
		Tree(std::string tstring,bool readBrLens=true);
		~Tree();
		
		Node *GetRoot() {return root; }
		unsigned GetNTax()	{return ntax; }
		unsigned GetMaxTaxonLabelLength(bool terminalsOnly=true);
		void SetName(std::string);
		std::string GetName();
		bool AttributesAreDirty()	{return attributesDirty; }
		void SetAttributesAreDirty(bool t) {attributesDirty=t; }
		Node **GetRecursiveNodeList();
		Node **GetRecursiveInternalNodeList();
		Node *GetTerminalNode(unsigned j) {
			assert(j < termTax.size());	
			return termTax[j];
		}
		void SetSimAttributes(SetOfTreeSimAttr *s) {
			lastUsedSimAtt = NULL;
			delete setOfSimAtt;
			setOfSimAtt=s;
		}
		
		SetOfTreeSimAttr &GetSimAttributesRef()	{
			assert (setOfSimAtt);
			return *setOfSimAtt;
		}
		SetOfTreeSimAttr * GetSimAttributesPtr() {
			return setOfSimAtt;
		}
		
		double LScoreAlreadyInitialized();
		double LScoreAlreadyInitialized(unsigned partnum,bool OnlyScoreDirty);
		void SimulateData(unsigned nRepsToDo);
		bool IsGood()
		{
			if(root)
				return root->IsGood();
			return false;
		}
		void CalcLikeBelow(Node *nod,Node **localRecNodList);
		void LikelihoodSweepOverAPartitionModel(Node **);
		double LScoreDirtyParts();
		void NotifyModelOfChangedParameter(PartModIndex pmi,Parameter *p);		
		void NotifyModelOfChangedParameter(PartModIndex pmi,FreqParamGroup *f);
		void NotifyModelsThatFreqParamChangesAreDone(PartModIndex pmi,FreqParamGroup *f);
		void PrepareLikeAttrAndAllNodes(unsigned partnum,unsigned modn);
		void PrepareLikeAttrForCalcuation(unsigned partnum,unsigned modn,Node *tempno);
		void PrepareSimAttrAndAllNodes(unsigned partnum,unsigned modn);
		double BranchLengthSmoothingPass(PartModIndex owner, std::vector<PartModIndex>* affected,double maxPasses,double delta,unsigned &NumBranchesAlreadyOptimized);
		void OneBranchLinearBrent(Node *,PositiveParameter *blen,double maxPasses,double delta);
		void OneBranchLinearBrent(double *condLikesBelow,Node *,double maxPasses,double delta);
		void BracketBranchLen(Node *nod,PositiveParameter *blen,double *a,double *ascore,double *b,double *bscore,double *c,double *cscore);
		void GetBracketNearMin(Node *nod,PositiveParameter *blen,double *a,double *ascore,double *b,double *bscore,double *c,double *cscore);
		void GetBracketNearMax(Node *nod,PositiveParameter *blen,double *a,double *ascore,double *b,double *bscore,double *c,double *cscore);
		void GetSecondPoint(Node *nod,PositiveParameter *blen,double *firPt,double *secPt,double *firSc,double *secSc,bool OnlyDecrease);
		double ScoreAffectedBranch(Node *,PositiveParameter *blen);
		void MultiplyBranches(double y);
		void Print(std::ostream &usrstream,bool withBrLen,bool tabbing=false,bool withReturns=true);
		void PrintBranchesWithMuts(std::ostream &usrstream, bool withBrLen,bool nonSynMut);
		void PrintModel(std::ostream &usrstream);
		void PrintTree(unsigned maxWidth,unsigned maxHeight,unsigned displayMode/*=branchDisplay(allNames)*/,unsigned optArg/*=0*/);
		void PrepareParsAttrForCalcuation(unsigned partnum);
		unsigned FirstPass(unsigned charact);
		void SecondPass(unsigned charact);
		double GetPStateFreq(unsigned charact,unsigned state,unsigned mode); //0 is maximize the character state, 1 is minimize it
		void NameInternalNodes(void);
		unsigned GetDepthTimesNSib(void) {return root->GetDepthTimesNSib(); }
		unsigned GetNumBranches();
		void SetArrayToAboveBottOfNode(double *cl,Node *tempno,unsigned *nStatesInThisSite);
		void ReRoot(Node *futldes);
		void RecodeTerminalsAndLikeAttr(unsigned p,unsigned m);
		static Tree *CreateYuleTree(unsigned ntaxInYule); //allocates
		void PrintTaxaBlock(std::ostream & outp);
		double GuessNumImpossibleBranches(unsigned slen) {assert (slen>100);	 double x=0.0;	root->GuessNumImpossibleBranches(slen,x);	return x; }



	protected:
		
		Node **recursiveNodeList;
		Node **recursInternalNodeList;
		unsigned ntax, nchar;
		std::vector<Node *> termTax;
		std::string name;
		unsigned memoryAllocation;
		Node *root;
		
		unsigned numTempArrsNeeded;
	
		bool attributesDirty;
		
		//TreeLikeAttr *lastScoredLikeAtt;
		//SetOfTreeLikeAttr *setOfLikeAtt;
		//TreeParsAttr *lastScoredParsAtt;
		//SetOfTreeParsAttr *setOfParsAtt;
		TreeSimAttr *lastUsedSimAtt;
		SetOfTreeSimAttr *setOfSimAtt;
		
		
};

} // namespace bull

#endif
