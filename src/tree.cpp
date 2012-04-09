
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

#include "tree.hpp"

#include <algorithm>

#include "ncl/nxstaxablock.h"
#include "string_extensions.hpp"

#define NUMBERNODES
using namespace std;
using namespace bull;

void Tree::PrintBranchesWithMuts(
  std::ostream &usrstream,
  bool withBrLen,
  bool nonSynMut) 
{
	usrstream<<"TREE "<< NxsString::GetEscaped(name) << " = [&U] ";
	if (withBrLen)
		usrstream << setprecision(8);
	if (root)
		root->PrintBranchesWithMuts(usrstream, withBrLen, nonSynMut);
	usrstream<<";";
}

void Tree::NameInternalNodes()
{	assert(root);
	assert(ntax);
	root->NameSelfAndInternalDescendants(ntax);
}


std::string Tree::GetName()
{	return name;
}

void Tree::PrepareSimAttrAndAllNodes(unsigned partnum,unsigned modn)
{
	setOfSimAtt->InitializeLikeAttrStatics(partnum,modn);
	for (unsigned i = 0; recursiveNodeList[i] != 0L; ++i)
		recursiveNodeList[i]->SetSimAtt(partnum,modn);
	lastUsedSimAtt = setOfSimAtt->GetSimAtt(partnum,modn);
	lastUsedSimAtt->GetModel()->AlertSharedMemory();
}


void Tree::SetName(std::string n)
{	ToUpper(n);
	name=n;
}
Tree::Tree(std::string tstring,bool readBrLens/*=true*/){
	double y;
	ntax=0;
	nchar=0;
	std::string::iterator c=tstring.begin();
	assert(*c == '(');
	root=new Node;
	Node *tempnode=root;
	hasBranchLengths=false;
	while (c!=tstring.end()){
		std::string tk;
		tk=GetTreeToken(tstring,c);
		/*	Paul's NCl trims the semicolon off the end
		if (c == tstring.end())
			{if(tk!=";")	throw TreeEncodingException("Expecting semicolon");
			}
		else 
		*/
		if (tk == "(") tempnode=tempnode->MakeDes();
		else if (tk == ")") tempnode=tempnode->GetAnc();
		else if (tk ==  ",")
			{
			tempnode=tempnode->GetAnc();
			tempnode=tempnode->MakeDes();
			}
		else if (tk == ":"){
			tk=GetTreeToken(tstring,c);
			if (readBrLens)
				{y=atof(tk.c_str());
				hasBranchLengths=true;
				tempnode->SetBranchLengthFromFile(y);
				}
			}
		else{
			UnderscoresToBlanks(tk);
			tempnode->SetName(tk);
			termTax.push_back(tempnode);
			ntax++;
			}
		}
	if (hasBranchLengths)
		{
		if (!root->NoIllegalBranchLengths())
			{cout<<"Setting all undefined or negative branches to 0.0"<<endl;
			root->SetAllIllegalBranches(0.0);
			hasBranchLengths=false;
			}
		root->CreateSetOfLikeAttrWithInputBrLensRecurs();
		}
	
	attributesDirty=true;
	recursiveNodeList=NULL;
	recursInternalNodeList=NULL;
	setOfSimAtt=NULL;
	tempOneBrentScored=0;
#ifdef NUMBERNODES
	//TEMPORARY
	NameInternalNodes();
#endif
	
	}



Tree::~Tree()
{	termTax.erase(termTax.begin(),termTax.end());
	delete root;
	delete [] recursiveNodeList;
	delete [] recursInternalNodeList;
	delete setOfSimAtt;
	
}
		

void Tree::PrintTaxaBlock(ostream &outp)
{	outp<<"begin taxa ; \ndimensions ntax="<<ntax<<"; \ntaxlabels ";
	for (unsigned j=0; j < ntax; j++)
		outp<<termTax[j]->name<<" ";
	outp<<"; \n";
}

void Tree::Print(ostream &usrstream,bool withBrLen,bool tabbing/*=false*/,bool withReturns/*=true*/)
{	usrstream<<"TREE "<<BlanksToUnderscores(name)<<"= [&U] ";
	if (tabbing && withReturns) usrstream<<"\n";
	if (withBrLen)
		usrstream<<setprecision(8);
	if (root)
		root->Print(usrstream,withBrLen,0,tabbing);
	usrstream<<"; ";
	if (withReturns)
		usrstream<<endl;
}
Node **Tree::GetRecursiveNodeList()
{	assert(root && ntax>1);
	if (!recursiveNodeList)
		{recursiveNodeList=new Node *[2*ntax+1];
		for (unsigned i=0; i < 2*ntax; i++)
			recursiveNodeList[i]=NULL;
		root->FillRecursiveNodeList(recursiveNodeList);
		}
	if (!recursInternalNodeList)
		{recursInternalNodeList=new Node *[ntax+1];
		for (unsigned i=0; i < ntax; i++)
			recursInternalNodeList[i]=NULL;
		root->FillRecursiveInternalNodeList(recursInternalNodeList);
		}
	return recursiveNodeList;
}


//should probably be done as SetOfTreeSimAttr function
void Tree::SimulateData(unsigned nRepsToDo)
{	
	Node **tempRNL = GetRecursiveNodeList();
	unsigned nd_num = 0;
	while (tempRNL[nd_num])
		++nd_num;
	assert(nd_num > 0);
	const unsigned root_nd_num = nd_num - 1;
	for (unsigned partn = 0; partn < setOfSimAtt->GetNParts(); partn++) {
		nd_num = root_nd_num;
		Node *ndp = tempRNL[nd_num];
		assert(ndp);	
		assert(ndp == root);
		unsigned numSimChars=nRepsToDo*setOfSimAtt->nChar[partn];
		assert(setOfSimAtt->GetNModels(partn) == 1); //other wise we'll have to pick a model for each character
		PrepareSimAttrAndAllNodes(partn,0);
		
		//generate the root sequence
		TreeSimAttr *tsa=setOfSimAtt->GetSimAtt(partn,0);
		
		SimNodeLikeAttr * siminf = ndp->simInfo;
		assert(siminf);
		assert(siminf->charsAsShorts.size() >= numSimChars);
		tsa->GenerateRandomSequence(numSimChars, &(siminf->charsAsShorts[0]));
		if (tsa->HasRateHet())
			tsa->GenerateNewSetOfRates(numSimChars);
		
		//sweep up the tree by backing up through the recursiveNodeList
		while (nd_num > 0) {
			nd_num--;
			ndp = tempRNL[nd_num];
			assert(ndp);
			siminf = ndp->simInfo;
			assert(siminf);
			assert(siminf->charsAsShorts.size() >= numSimChars);
			short *tnC = &(siminf->charsAsShorts[0]);

			assert(ndp->anc);
			SimNodeLikeAttr * ancsiminf = ndp->anc->simInfo;
			assert(ancsiminf);
			assert(ancsiminf->charsAsShorts.size() >= numSimChars);
			const short *anC = &(ancsiminf->charsAsShorts[0]);
			
			if (tsa->HasRateHet())
				siminf->SimulateSeq(numSimChars, anC, tsa->GetRatesPtr());
			else
				siminf->SimulateSeq(numSimChars, anC);
			
			if (nRepsToDo == 1 && ndp->ldes) {
				if (partn == 0) {
					ndp->hasNucMutation = false;
					ndp->hasNonSynMutation = false;
				}
				for (unsigned ij=0; ij < setOfSimAtt->nChar[partn]; ij++) {
					if (tnC[ij] != anC[ij]) {
#						warning "dangerous code assumes SSRFCodonSubMod model here-temporary"
						ndp->hasNucMutation = true;
						SSRFCodonSubMod * ssrfmod = dynamic_cast<SSRFCodonSubMod *>(siminf->GetModel());
						assert(ssrfmod);
						unsigned tnAA = ssrfmod->codLocToAminoAcidLoc.at(tnC[ij]);
						unsigned ancAA = ssrfmod->codLocToAminoAcidLoc.at(anC[ij]);
						if (tnAA != ancAA) {
							ndp->hasNonSynMutation = true;
							break;
						}
					}
				}
			}
		}
	}
}


unsigned Tree::GetMaxTaxonLabelLength(bool terminalsOnly/*=true; */)
{	
	unsigned maxlen = 0;
	if (terminalsOnly) {
		for (vector < Node *>::iterator ttIt= termTax.begin(); ttIt!=termTax.end(); ttIt++)
			if ((*ttIt)->name.length()>maxlen)
				maxlen=(*ttIt)->name.length();
	}
	else {
		Node **tempnl=GetRecursiveNodeList();
		while (*tempnl) {
			if ((*tempnl)->name.length()>maxlen)
				maxlen=(*tempnl)->name.length();
			tempnl++;
		}
	}
	return maxlen;
}
	



