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
 

#ifndef _NODE_
#define _NODE_


#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 

#include "encoded_chars.hpp"
#include "like_attribute_sets.hpp"
#include "util.hpp"

#define BLOSPSEUDTERM 4
#define BLOSFINISHED 1
#define BLOSBELOWSUBTREE 2


namespace bull {

class NoFreeWorkArrayException {};
class RequestedWorkArrayNotFoundException {};
class Node;

void MakeFirstArrayProduct(double *p,double *f,double *s,int lenOfArray);
void SetFirstArrayEqualToSec(double *f,double *s,int lenOfArray);

class Node{
	friend class Tree;
	//friend class BullShell;
	private :
		std::string name;
		double *blen;
		SetOfLikeAttr *setOfLikeInfo;
		//SetOfParsAttr *setOfParsInfo;
		Node *ldes;
		Node *rdes;
		Node *next;
		Node *anc;
	
	public:

	bool hasNucMutation; //temporary simulation hack to identify nodes that have a true synapormorphy  ONLY WORKS if Simulating ONE DATA SET AT A TIME
	bool hasNonSynMutation; //temporary simulation hack to identify nodes that have a true synapormorphy  ONLY WORKS if Simulating ONE DATA SET AT A TIME

	SimNodeLikeAttr *simInfo;
	void ReplaceSetOfLikeInfoButKeepBranchLen(SetOfLikeAttr *LI);
	void SetSetOfLikeAttr( SetOfLikeAttr *swapSOLA) { setOfLikeInfo=swapSOLA; }
	SetOfLikeAttr *GetSetOfLikeAttr( void ) {return	 setOfLikeInfo; }
	Node();
	Node(bool withBranchLenField);
	void SharedNodeConstruction();
	~Node();
	Node *CopyToShare();
	Node *Copy();
	bool IsTerm()	{return !ldes; }
	Node *GetAnc();
	int GetNchar();
	void WriteTermTax(std::vector<Node *> *);
	void WriteIfTerm(std::vector<Node *> *);
	Node *GetRoot();
	int GetNtax();
	int GetNtax(int);
	void AddDes(Node *d);
	Node *MakeDes();
	Node *MakeDes(bool withBLField);
	void SetBranchLengthFromFile(double); //used in reading the branches from a file (they are stored in the node's blen temporarily as opposed to a like attr
	void SetName(std::string);
	std::string GetName();
	void SetAllIllegalBranches(double len);
	bool NoIllegalBranchLengths();
	void CalcVariableCondLikelihoodTopDown();
	void CondLikeArrayMultPrChanges(double *cl);
	double GetLogLikelihood();
	Node **FillRecursiveNodeList(Node **);
	Node **FillRecursiveInternalNodeList(Node **rnl);
	void SetSimAtt(int part,int mod)	{ 
		assert(setOfLikeInfo);
		simInfo=(SimNodeLikeAttr *)setOfLikeInfo->GetLikeAtt(part,mod); 
		//TEMP
		blen=simInfo->GetBLenPtr();
		}
	void CreateSetOfLikeAttrWithInputBrLensRecurs();
	bool IsGood(void);
	Node *GetNextNodeToOptBranchOn();
	void MultiplyBranch(double y);
	void MultiplyBranchRecurs(double y);
	void AddToBranchLengthRecurs(double y);
	void AddToBranchLength(double y);
	void PrintBranchesWithMuts(std::ostream &usrstream,bool withBrLen, bool nonSyn);
	void Print(std::ostream &usrstream,bool withBrLen,int level=0,bool tabbing=false);
	int NameSelfAndInternalDescendants(int n);
	void ShowNode(int windowWidth, std::vector<int> prevVert);
	int GetNumDes(int n,bool termOnly);
	int MaxNumInternalNodesToTip();
	bool MustBeAChangeOnBranchUppass(int);
	int GetDepthTimesNSib();
	void DirectToNewAnc(Node* newAnc);
	void SwapBranchLengths(Node *sourc);
	void AddBranchLengths(Node *sourc);
	int GetNumBranchesAbove();
	int GetHeightAbove(int height);
	void EvolveRateAndChangeBranchLengths(double rateAtBeg,double roeoroe,double *sumBrLens,double *sumOfBrLenTimes);
	void EvolveRateAndChangeBranchLengths(double rateAtBeg,double roeoroe,double maxRate,double minRate,double *sumBrLens,double *sumOfBrLenTimes);

	void ResetLWAs();
	void FreeCondLikeAboveBottom();
	static int currPrioritySetting;
	static int NumBranchesOptimized;
	static int branchLengthOptLogNum;
	static double TreesCurrentLikelihood;

	int UpdateBranchesInRecursiveOrder(double *optBranchLengths, int nBranchesToUpdate);
	void SetMinimumTerminalBranchLength(double mb);
	double GetMeanBrLensOfDescendantClades(); //UltraDev
	void AssignHeight();
	void SetBranchLengthsAccordingToNodeHeights(double rootAge);
	void GuessNumImpossibleBranches(int slen,double &x);
	//void PrintBranchesInRecursiveOrder();
};

	double BranchLengthLinearBrent(double *branchLength,double *below,double *above,double *work,Node *currNode,int *nStatesInThisSite,Model **modArrays,double maxPasses,double delta,int numSitesInProtein);
	double ScoreBranchLength(double origBranch,double *below,double *above,double *work,Node *currNode,int *nStatesInThisSite,Model **modArrays,int numSitesInProtein);
	double ScoreZeroBranchLength(double *below,double *above,double *work,Node *currNode,int *nStatesInThisSite,Model **modArrays,int numSitesInProtein);
	void MakeBelowTopFromBelowBottom(double bLen,double *top,double *bottom,int *nStatesInThisSite,Model **modArrays,int numSitesInProtein);
	void GetBranchLengthBracket(int numPtsToTrust,double *a,double *b, double *c, double *ascore, double *bscore,double *cscore,double *below,double *above,double *work,Node *currNode,int *nStatesInThisSite,Model **modArrays,int numSitesInProtein);
	void GetMiddleBranchLength(double *a,double *b, double *c, double *ascore, double *bscore,double *cscore,double *below,double *above,double *work,Node *currNode,int *nStatesInThisSite,Model **modArrays,int numSitesInProtein);
	void WriteBrLenLike(double *below,double *above,int *nStatesInThisSite,Model **modArrays,int numSitesInProtein,int nOfBranch);

} // namespace bull 
#endif
