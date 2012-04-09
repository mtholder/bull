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
 

#include "node.hpp"
#include "string_extensions.hpp"

using namespace bull;


std::string Node::GetName(){
	return name;
	}
Node::~Node(){
	//Node isn't responsible for deleting blen, the LikeAttribute will do that even if the node allocates the memory (it will be put into an LikeAttr by the Tree)
	delete ldes;
	delete next;
	delete setOfLikeInfo;
	}
void Node::SetAllIllegalBranches(double len)
{	assert(len >= 0.0);
	if (ldes)	ldes->SetAllIllegalBranches(len);
	if (next)	next->SetAllIllegalBranches(len);
	if (anc)
		{if(!blen)	SetBranchLengthFromFile(len);
		if (*blen < 0.0) *blen=len;
		}
}
Node *Node::MakeDes(){
	if (!ldes){
		ldes=new Node;
		rdes=ldes;
		}
	else{
		rdes->next=new Node;
		rdes=rdes->next;
		}
	rdes->anc=this;
	return rdes;
	}


Node *Node::GetAnc(){
	return anc;
	}
	

void Node::SetBranchLengthFromFile(double y){
	blen=new double;
	*blen=y;
	}


void Node::SetName(std::string n){
	name=n;
	}
	
bool Node::NoIllegalBranchLengths()
{	if (anc)
		{if(!blen)
			return false;
		if (*blen < 0.0)
			return false;
		}
	if (ldes)	
		if (!ldes->NoIllegalBranchLengths())
			return false;
	if (next)	
		if (!next->NoIllegalBranchLengths())
			return false;
	return true;
}	


int Node::NameSelfAndInternalDescendants(int n)
{	if (!ldes) return n;
	Node *tempno=ldes;
	while (tempno)
		{n=tempno->NameSelfAndInternalDescendants(n);
		tempno=tempno->next;
		}
	name="";
	name+=n;
	return n+1;
}
void Node::CreateSetOfLikeAttrWithInputBrLensRecurs()
{	if (anc) 
		assert(blen);
	assert(!setOfLikeInfo); //this should only be called to read branch lengths in from a file
	setOfLikeInfo= new SetOfLikeAttr(1);
	setOfLikeInfo->SetNumModsInPart(0,1);
	LikeAttr *templa;
	templa=new LikeAttr();
	if (blen)
		{templa->SetBLen(*blen);
		delete blen;
		blen=NULL;
		}
	else	templa->SetBLen(0.0);
	blen=templa->GetBLenPtr();
	setOfLikeInfo->AddLikeAtt(0,0,templa);
	if (ldes)	ldes->CreateSetOfLikeAttrWithInputBrLensRecurs();
	if (next)	next->CreateSetOfLikeAttrWithInputBrLensRecurs();
}



bool Node::IsGood(void)
{	
	if (ldes)	{if(!ldes->IsGood())	
					return false;
				if (!rdes)	
					return false;
				if (rdes == ldes) 
					return false;
				}
	if (anc && !next )
		if (this!=anc->rdes)
			return false;
	if (next)	if (!next->IsGood())	
					return false;
	if (rdes&&!ldes)	
		return false;
	return true;

}
	


//@ modified 6/7/02 after not looking at this code for a _really_ long time.
//	the assert was tripping as I destroyed simulation attributes (because I had simulated under a model with lots parts)
//	I stole code from the CreateSetOfLikeAttrWithBrLenFromFile function for the case when the argument is NULL
//

void Node::ReplaceSetOfLikeInfoButKeepBranchLen(SetOfLikeAttr *LI)
{	//should be called only for getting user input brlens, if modifying beware of deletion of LikeAttr which can cause double deletion of blen (if blen is shared between models)
	assert(setOfLikeInfo);
	if (LI)
		{
		assert(setOfLikeInfo->GetNParts() == 1 && setOfLikeInfo->GetNModels(0) == 1);
		for (unsigned p=0; p < LI->GetNParts(); p++)
			for (unsigned m=0; m < LI->GetNModels(p); m++)
				LI->GetLikeAtt(p,m)->SetBLenPtr(setOfLikeInfo->GetLikeAtt(0,0)->GetBLenParameterPtr());
		}
	else
		{
		LI = new SetOfLikeAttr(1);
		LI->SetNumModsInPart(0,1);
		LikeAttr *templa;
		templa=new LikeAttr();
		templa->SetBLenPtr(setOfLikeInfo->GetLikeAtt(0,0)->GetBLenParameterPtr());
		blen=templa->GetBLenPtr();
		LI->AddLikeAtt(0,0,templa);
		}
	setOfLikeInfo->GetLikeAtt(0,0)->DetachBLenPtr(); //so we don't destroy the branchlengths
	delete setOfLikeInfo; 
	setOfLikeInfo=LI;
}

Node **Node::FillRecursiveNodeList(Node **rnl)
{	if (ldes)	{rnl=ldes->FillRecursiveNodeList(rnl);
				Node *tempno=ldes->next;
				while (tempno)
					{rnl=tempno->FillRecursiveNodeList(rnl);
					tempno=tempno->next;
					}
				}
	*rnl++=this;
	*rnl=NULL;
	return rnl;
}

Node **Node::FillRecursiveInternalNodeList(Node **rnl)
{	
	if (ldes)	{rnl=ldes->FillRecursiveInternalNodeList(rnl);
				Node *tempno=ldes->next;
				while (tempno)
					{rnl=tempno->FillRecursiveInternalNodeList(rnl);
					tempno=tempno->next;
					}
				*rnl++=this;
				*rnl=NULL;
				}
	return rnl;
}



void Node::Print(
  std::ostream &usrstream,
  bool withBrLen,
  int level/*=0*/,
  bool tabbing/*=false*/) {
	if (tabbing)
		{/*if(!anc )	{usrstream<<"("; //add opening parentheses for the whole tree
					level++;
					}*/
		if (ldes)	{
					usrstream<<"\n";
					for (int i=0; i < level; i++) usrstream<<"\t";
					usrstream<<"(";
					level++;
					ldes->Print(usrstream,withBrLen,level,true);
					assert(ldes!=rdes);
					if (withBrLen && anc && blen)
						usrstream<<":"<<*blen;
					level--;
					}
		else		{std::string tempn=name;
					BlanksToUnderscores(tempn);
					usrstream<<"\n";
					for (int i=0; i < level; i++) usrstream<<"\t";
					usrstream<<tempn;
					if (withBrLen && blen)
						usrstream<<":"<<*blen;
					}
		if (next)	{usrstream<<",";
					next->Print(usrstream,withBrLen,level,true);	
					}
		else		{level--;
					if (anc)
						{usrstream<<"\n";
						for (int i=0; i < level; i++) usrstream<<"\t";
						usrstream<<")";
						}
					}
		}
	else
		{//if(!anc) usrstream<<"("; //add opening parentheses for the whole tree
		if (ldes)	{usrstream<<"(";
					ldes->Print(usrstream,withBrLen);
					assert(ldes!=rdes);
					if (withBrLen && anc && blen)
						usrstream<<":"<<*blen;
					}
		else		{std::string tempn=name;
					BlanksToUnderscores(tempn);
					usrstream<<tempn;
					if (withBrLen && blen)
						usrstream<<":"<<*blen;
					}
		if (next)	{usrstream<<",";
					next->Print(usrstream,withBrLen);	
					}
		else		{if(anc)	usrstream<<")";
					}
		}
}


void Node::PrintBranchesWithMuts(std::ostream &usrstream, bool withBrLen, bool nonSynMut)
{	
	Node *currChild = ldes;
	if (currChild) {
		const bool mutThatCounts = hasNucMutation && ((!nonSynMut) || hasNonSynMutation);
		if(mutThatCounts || anc == NULL)
			usrstream << '(';
		while (currChild) {
			currChild->PrintBranchesWithMuts(usrstream, withBrLen, nonSynMut);
			if (currChild != rdes)
				usrstream << ',';
			currChild = currChild->next;
		}
		if(mutThatCounts || anc == NULL) {
			usrstream << ')';
			if (withBrLen && anc && blen)
				usrstream << ':' << *blen;
		}
	}
	else {
		usrstream << NxsString::GetEscaped(name);
		if (withBrLen && blen)
			usrstream<<":"<<*blen;
	}
}


void Node::SharedNodeConstruction(){
	ldes=rdes=next=anc=0L;
	blen=NULL;
	setOfLikeInfo=NULL;
	simInfo=NULL;
}

Node::Node() {
	SharedNodeConstruction();
	blen=NULL;
}



