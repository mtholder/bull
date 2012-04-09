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
 
#include<iostream>

#include "model.hpp"
#include "lin_alg.hpp"
#include "tools.hpp"
#include "matrices.hpp"
#include "util.hpp"

using namespace bull;
using std::cout;

void Model::PrintPAUPLsetCommand() const {
	cout << "lset pinv = ";
	if (pinv)
		cout << pinv->val;
	else
		cout<<"0.0";
	cout << " rates = ";
	if (ngamcat > 1)
		cout << "gamma shape = " << gammashape->val << " ";
	else
		cout << "equal ";
	}

double RateManager::defPinv(.1);

double RateManager::defGammaShape(.5);

ModelWEig::~ModelWEig()
{
	delete [] REigenValues;
	delete [] EigInvEigMat;
	delete [] ImEigenValues;
	delete [] EigExp;
	free_psdmatrix(REigenVector);
	free_psdmatrix(InvEigenVector);
	free_pscmatrix(CEigenVector);
	free_pscmatrix(CInvEigenVector);
	free_psdmatrix(qMatrix);
}

void ModelWEig::SharedConstruction(unsigned n)
{
	REigenValues = new double[n];
	ImEigenValues = new double[n];
	REigenVector = psdmatrix(n);
	InvEigenVector = psdmatrix(n);
	CEigenVector = pscmatrix(n);
	CInvEigenVector = pscmatrix(n);
	qMatrix = psdmatrix(n);
	EigExp = new double[n];
	eigencalc = false;

}




void ModelWEig::UpdatePMatrix(double **pmats,double blen)	
{	
	assert(REigenValues && ImEigenValues && REigenVector && InvEigenVector && CEigenVector && CInvEigenVector);
	if (!eigencalc) {
		CalculateQ();
		IsComplex = GetEigens(nstates,qMatrix,REigenValues,ImEigenValues,REigenVector,InvEigenVector,CEigenVector,CInvEigenVector);
#		ifdef	CONDENSE_MATRICES
			CalculateAndCondenseEigInvEigMult(*REigenVector,InvEigenVector,nstates,EigInvEigMat);
#		else
#			error "linalg assumes that CONDENSE_MATRICES is 1
			CalculateGlobalEigInvEigMult(*REigenVector,InvEigenVector,nstates,EigInvEigMat);
#		endif
		eigencalc = true;
	}
	if (!IsComplex)
		ChangeMatrix(blen,*pmats,nstates,REigenValues,EigInvEigMat,EigExp);
	else if(ComplexChangeMatrix(nstates,pmats,blen,REigenValues,ImEigenValues,CEigenVector,CInvEigenVector))
		throw LinAlgProb("Error in ComplexChangeMatrix");
}

//this is a faster version of UpdatePMatrix that you can use if only one column of the PMatrix is needed
void ModelWEig::UpdatePMatrix(double **pmats,double blen,unsigned onlycolumn)	
{
	assert(REigenValues && ImEigenValues && REigenVector && InvEigenVector && CEigenVector && CInvEigenVector);
	if (!eigencalc) {
		CalculateQ();
		IsComplex = GetEigens(nstates,qMatrix,REigenValues,ImEigenValues,REigenVector,InvEigenVector,CEigenVector,CInvEigenVector);
#		ifdef	CONDENSE_MATRICES
			CalculateAndCondenseEigInvEigMult(*REigenVector,InvEigenVector,nstates,EigInvEigMat);
#		else
			CalculateGlobalEigInvEigMult(*REigenVector,InvEigenVector,nstates,EigInvEigMat);
#		endif
		eigencalc = true;
	}
	if (!IsComplex) {
#	ifdef CONDENSE_MATRICES
		if (onlycolumn)
			ChangeColumn(blen,*pmats,nstates,REigenValues,EigInvEigMat,EigExp,onlycolumn);
		else	//can't get a speed up if the matrix is condensed and the only column is 0 because this the column obtained by subtraction
			ChangeMatrix(blen,*pmats,nstates,REigenValues,EigInvEigMat,EigExp);
#	else
		ChangeColumn(blen,*pmats,nstates,REigenValues,EigInvEigMat,EigExp,onlycolumn);
#	endif
	}
	else if(ComplexChangeMatrix(nstates,pmats,blen,REigenValues,ImEigenValues,CEigenVector,CInvEigenVector))
		throw LinAlgProb("Error in ComplexChangeMatrix");
}

//this is a faster version of UpdatePMatrix that you can use if only one column of the PMatrix is needed
void ModelWEig::UpdatePRow(double **pmats,double blen,unsigned onlycolumn)	
{	
	assert(REigenValues && ImEigenValues && REigenVector && InvEigenVector && CEigenVector && CInvEigenVector);
	if (!eigencalc) {
		CalculateQ();
		IsComplex = GetEigens(nstates,qMatrix,REigenValues,ImEigenValues,REigenVector,InvEigenVector,CEigenVector,CInvEigenVector);
#		ifdef	CONDENSE_MATRICES
			CalculateAndCondenseEigInvEigMult(*REigenVector,InvEigenVector,nstates,EigInvEigMat);
#		else
			CalculateGlobalEigInvEigMult(*REigenVector,InvEigenVector,nstates,EigInvEigMat);
#		endif
		eigencalc = true;
	}
	if (!IsComplex) {
#		ifdef CONDENSE_MATRICES
			ChangeRow(blen,*pmats,nstates,REigenValues,EigInvEigMat,EigExp,onlycolumn);
#		else				
			ChangeRow(blen,*pmats,nstates,REigenValues,EigInvEigMat,EigExp,onlycolumn);
#		endif		
	}
	else if(ComplexChangeMatrix(nstates,pmats,blen,REigenValues,ImEigenValues,CEigenVector,CInvEigenVector))
		throw LinAlgProb("Error in ComplexChangeMatrix");
}

ModelWEig::ModelWEig(unsigned n,double *preAllocEIEM)	
		: Model(n) 
{
	SharedConstruction(n);
	EigInvEigMat = preAllocEIEM;
}

//I can't remember where I pinched this code from
void RateManager::CalculateRates(void )
{	
	assert(gammashape && !rates.empty() && ngamcat>1);
	unsigned i;
	double a, b, invK, x, y;
	double N = (double)ngamcat;
	double twoShape = 2.0 * gammashape->val;

	invK = 1.0 / N;
	b = 0.0;
	x = 1.0;  /* (in case user sets nRateCategs=1) */
	for (i = 0; i < ngamcat - 1; i++) {
		a = b;
		b = PointChi2( invK * (i+1), twoShape ) / twoShape;
		x = gammq( gammashape->val + 1.0, b * gammashape->val );
		if ( a > 0.0 )
			y = gammq( gammashape->val + 1.0, a * gammashape->val );
		else
			y = 1.0;
		rates[i] = (y - x) * N;
	}
	rates[ngamcat - 1] = x * N;
}


//default no rate het
RateManager::RateManager() {
	ngamcat=1;
	pinv=NULL;
	gammashape=NULL;
}

void RateManager::InitializeParameters()
{	
	if (gammashape) {
		if(gammashape->StartWithCurrent())
			{
			}
		else if (gammashape->StartWithRandom())
			throw IncompleteModel("Random Function to initialize parameters isn't available yet");
		else if (gammashape->StartWithApproximation())
			throw IncompleteModel("Initial approximation of parameters isn't available yet");
		else if (gammashape->StartWithDefault())
			gammashape->SetToDefault();
		else
			throw BadSettings("No starting value of a parameter has been defined");
		}
	if (pinv) {
		if(pinv->StartWithCurrent())
			{}
		else if (pinv->StartWithRandom())
			throw IncompleteModel("Random Function to initialize parameters isn't available yet");
		else if (pinv->StartWithApproximation())
			throw IncompleteModel("Initial approximation of parameters isn't available yet");
		else if (pinv->StartWithDefault())
			pinv->SetToDefault();
		else
			throw BadSettings("No starting value of a parameter has been defined");
	}
}


Model::Model(unsigned n) 
	:RateManager(),
	pMatrix(NULL)
{
	nstates = n;
	pMatrix = psdmatrices(n,1);
	pmatcalc = false;
	stateFreqs = NULL;
}

Model::~Model(){
	free_psdmatrices(pMatrix, ngamcat);
	for (std::vector<Parameter *>::iterator pIt = ownedParamVec.begin(); pIt != ownedParamVec.end(); ++pIt) {
		Parameter *pptr = *pIt;
		delete pptr;
	}
	delete stateFreqs;
	
}

RateManager::~RateManager()
{	
	delete gammashape;
	delete pinv;
}



