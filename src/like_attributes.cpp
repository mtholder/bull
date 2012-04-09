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
 

#include "like_attributes.hpp"
#include "bull.hpp"


using namespace bull;


int LikeAttr::currNChar(1);
int LikeAttr::currNStates(4);
int LikeAttr::currShortPerChar(1);
int LikeAttr::currNRateCats(1);
int LikeAttr::currNStatesInLastShort(1);
bool LikeAttr::currModelDirty(true);
double LikeAttr::Multiplier(1.53);
bool SimNodeLikeAttr::verboseSimulation = false;

#if defined (CHAR_BY_CHAR) && CHAR_BY_CHAR
	int LikeAttr::currCharIndex(0);
#endif

LikeAttr::LikeAttr()//Constructor for internal node's attributes
{	model=NULL;
	Pmat=NULL;
	blen=NULL;
	blenMod=NULL;
}

void LikeAttr::SetBLen(double d)
{	if (blen)
		blen->SetCurrent(d);
	else {
		const int r = par(CUR)|par(MIN);
		blen= new BoundedParameter(d, 0.0, 0.0, r, bull::DEFAULT_START_EDGE_LEN);
	}
}	


LikeAttr *LikeAttr::Copy()
{	LikeAttr *tempLA;
	if (model)
		tempLA=new LikeAttr(model);
	else
		tempLA=new LikeAttr();
	if (blen)
		tempLA->blen=new BoundedParameter(blen->val,0.0,0.0,blen->GetSetting(),
											bull::DEFAULT_START_EDGE_LEN);
	//TEMPORARY COPY Leaves blen mod uncopied  MEMORY LEAK old blen still allocated I don't know why anymore!!!
	return tempLA;
}



TreeSimAttr::TreeSimAttr(
	unsigned nc,
	const Model &m)
	:model(m),
	rates(nc, 1.0)
{
	ownsBrLens = false;

}


Model* LikeAttr::GetModel() {
	return model;
}
LikeAttr::~LikeAttr()
{	//assert(SharePmatMemory());
	delete blen;
	delete blenMod;
}

void TreeSimAttr::GenerateRandomSequence(unsigned nc, short *dest) const
{	
	assert(nc <= GetNChar());
	const double * const * freq = GetStateFreqs() ;
	const unsigned ns = GetNStates();
	for (unsigned i = 0; i < nc; i++) {
		const short state_ind = (short) GetRandomIndexFromFreqs(ns, freq);
		dest[i] = state_ind;
	}
}

void TreeSimAttr::GenerateNewSetOfRates(unsigned nc)
{	
	const double pnv=model.GetPInv();
	rates.resize(nc);
	if (pnv>0.0) {
		if(model.GetNRateCats() > 1) {
			double sp=model.GetShapeParam();
			for (unsigned i=0; i < nc; i++)
				rates[i] = (RandomNumber() < pnv ? 0.0 : (rndgamma(sp))/sp);
		}
		else {
			for (unsigned i=0; i < nc; i++)
				rates[i] = (RandomNumber() < pnv ? 0.0 : 1.0);
		}
	}
	else {
		assert(model.GetNRateCats()>1);
		double sp=model.GetShapeParam();
		for (unsigned i=0; i < nc; i++)
			rates[i] = (rndgamma(sp))/sp;
	}
}

LikeAttr::LikeAttr(Model *m)//Constructor for internal node's attributes
{	model=m;
	Pmat=m->GetPmat();
	blen=NULL;
	blenMod=NULL;
}


SimNodeLikeAttr::SimNodeLikeAttr (unsigned nc,Model *m)
	: LikeAttr(m),
	charsAsShorts(nc, 0)
{
}

void SimNodeLikeAttr::SimulateSeq(unsigned nc, const short *ancSeq, const double *rates)
{	//rate heterogeneity
	assert(model);
	assert(nc <= this->charsAsShorts.size());
	short* cis = &charsAsShorts[0];
	double transformedBranch=(blen->val)/(1.0-model->GetPInv());
	double ** tpmat = *Pmat;
	if (model->GetNRateCats()>1)
		{
		for (unsigned i=0; i < nc; i++)//rates could be anything
			{
			const short a = ancSeq[i];
			if (rates[i] > 0.0) {
				model->UpdatePRow(tpmat, (rates[i])*transformedBranch, a);
				cis[i] = (short) GetRandomIndexFromFreqs(currNStates,tpmat[a]);
			}
			else
				cis[i] = a;
		}
	}
	else {
		model->UpdatePMatrix(*Pmat, transformedBranch); //get whole matrix, the sites will either be rate 1 or 0
		for (unsigned i=0; i < nc; i++) {
			const short a = ancSeq[i];
			cis[i] = (rates[i] > 0.0 ? (short) GetRandomIndexFromFreqs(currNStates,tpmat[a]) : a);
		}
	}
			
}	

void SimNodeLikeAttr::SimulateSeq(unsigned nc, const short *ancSeq)
{	//no rate heterogeneity
	assert(model);
	assert(nc <= this->charsAsShorts.size());
	double ** tpmat = *Pmat;
	model->UpdatePMatrix(tpmat, blen->val/(1.0 - model->GetPInv()));
	
	if (SimNodeLikeAttr::verboseSimulation) {
		std::cout << "brlen = " << blen->val << " \n";
		std::cout << "-1";
		for (unsigned s = 0; s < currNStates; ++s)
			std::cout << ' ' << s ;
		std::cout << '\n';
		for (unsigned f = 0; f < currNStates; ++f) {
			std::cout << f;
			for (unsigned t = 0; t < currNStates; ++t) {
				std::cout << ' ' <<  tpmat[f][t];
			}
			std::cout << '\n';
		}
		std::cout.flush();
	}
	
	for (unsigned i=0; i < nc; i++) {
		const short a = ancSeq[i];
		const short d = (short) GetRandomIndexFromFreqs(currNStates, tpmat[a]);
		charsAsShorts[i] = d;
	}
}





