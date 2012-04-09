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
 

#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>

#include "basic_bull.hpp"

#include <map>
#include "encoded_chars.hpp"
#include "char_encoding.hpp"
#include "mth_exception.hpp"
#include "parameter.hpp"
#include "complex.hpp"

namespace bull {

class BadParams: public MTHException 
{
	public:
		BadParams(const char *p) 
			:MTHException(p) 
			{
			}
};

class LinAlgProb: public MTHException 
{
	public:
		LinAlgProb(const char *p) 
			:MTHException(p) 
			{
			}
};


class IncompleteModel: public MTHException 
{
	public:
		IncompleteModel(const char *p) 
			:MTHException(p) 
			{
			}
};

class BadSettings: public MTHException 
{
	public:
		BadSettings(const char *p)
			:MTHException(p) 
			{}
};


class RateManager
{
	public :
		static void set_default(double p, double g) {
			defPinv = p;
			defGammaShape = g;
		}

		RateManager(); //default no rate het
		RateManager(unsigned, double); //gamma rates with shape
		RateManager(double); //pinv
		RateManager(double, unsigned, double); //pinv + gamma shape not specified
		
		virtual ~RateManager();
		virtual void InitializeParameters();
		
		unsigned GetNRateCats() const {
			return ngamcat;
		}
		
		double GetPInv() const {
			if (pinv)
				return pinv->val;
			return 0.0;
		}
		
		virtual unsigned GetNParams() const {
			unsigned i = (this->pinv ? 1 : 0);
			if (this->ngamcat > 1)
				++i;
			return i;
		}	
		
		virtual Parameter *GetParameter(unsigned i) {
			assert(i == 0 ||i == 1);
			assert(this->pinv || this->gammashape);
			if (i == 0)
				return (this->pinv ? this->pinv : this->gammashape);
			assert(this->pinv); //if i=1 there must be both pinv and gammashape
			return this->gammashape;
		}

		bool HasRateHet() const {
			if (this->ngamcat > 1)
				return true;
			return (this->pinv && (this->pinv->val > 0.0));
			}
		
		double GetShapeParam() const {
			assert(this->gammashape);
			return gammashape->val;
		}

	protected:
		void CalculateRates();

		const static unsigned maxcat=20;
		static double defPinv;
		static double defGammaShape;
		
		std::vector<double> rates;
		unsigned ngamcat;
		Parameter *pinv;
		Parameter *gammashape;
};

class Model : public RateManager
{
	public :
	
		Model(unsigned n);
		Model(unsigned n, double ***pMat); //for memory sharing (only SSRF now)
		Model(unsigned n, double pi); //some invariant sites
		Model(unsigned n, unsigned ncats, double gamsh); //gamma rates
		Model(unsigned n, double pi, unsigned ncats, double gamsh); //gamma +pinvariant
		virtual ~Model();

		virtual unsigned NumStatesInLastShort() {
			//written to allow odd coding of SSRFCodonSubModel which overrides
			throw	MTHException("Entered unwritten code NumStatesInLastShort" );
		}
		virtual unsigned NumShortsPerCharacter() {
			//written to allow odd coding of SSRFCodonSubModel which overrides
			throw	MTHException("Entered unwritten code NumShortsPerCharacter" );
		}
		virtual void EncodeACharacter(short *, short *, unsigned , bool ) {
			//written to allow odd coding of SSRFCodonSubModel which overrides
			throw	MTHException("Entered unwritten code Model::EncodeACharacter" );
		}

		virtual void AlertSharedMemory() const {
			//written to allow odd coding of SSRFCodonSubModel which overrides
			throw	MTHException("Entered unwritten code AlertSharedMemory" );
		}
		virtual double GetMultiplier()	const {
			throw	MTHException("Entered unwritten code GetMultiplier" );
		}
		virtual void SetMultiplier(double) {
			throw	MTHException("Entered unwritten code SetMultiplier" );
		}

		virtual void PrintPAUPLsetCommand() const;
		
		double ***GetPmat() {
			return pMatrix;
		}
		
		unsigned GetNParams() const {
			return nparams + RateManager::GetNParams();
		}
		
		virtual void UpdatePMatrix(double **pm, double b)=0;
		
		virtual void UpdatePMatrix(double **pm, double b, unsigned /*onlycol*/) {
			UpdatePMatrix(pm, b);
		}
		
		virtual void UpdatePRow(double **pm, double b, unsigned /*onlycol*/) {
			UpdatePMatrix(pm, b);
		}
		
		virtual void UpdatePmat(double b) {
			if (ngamcat>1) { 
				for (unsigned i=0; i < ngamcat; ++i)
					UpdatePMatrix(pMatrix[i], rates[i]*b);
			}
			else
				UpdatePMatrix(*pMatrix, b);
		}
		
		virtual void UpdatePmat(double b, unsigned onlycol)	{
			if (ngamcat>1) { 
				for (unsigned i=0; i < ngamcat; i++)
					UpdatePMatrix(pMatrix[i], rates[i]*b, onlycol);
			}
			else 
				UpdatePMatrix(*pMatrix, b, onlycol);
		}
		
		
		Parameter *GetParameter(unsigned i) {
			//should be overwritten if the model has multiple freqgroups
			if (i < nparams)
				return allParamVec[i];
			return RateManager::GetParameter(i - nparams);
		}
		

		virtual const double * const * GetStateFreqsConst() const {
			assert(stateFreqs);
			return stateFreqs->GetStateFreqsConst();
		}
	
		inline	unsigned GetNFreqGroups() const {
			return nfreqParamGroups;
		}
		
		virtual unsigned GetEncodingType() const = 0;
		
		virtual unsigned GetNStates() const = 0;
		
		virtual void InitializeParameters()=0;
	
		virtual void ParameterHasChanged(Parameter *p) {
			if (p == this->gammashape)
				CalculateRates();
		}
		
		virtual void ParameterHasChanged(FreqParamGroup *)	{
			throw MTHException("Model::ParameterHasChanged(FreqParamGroup *p) has been called");
		}
		
		virtual void FreqParamChangesShouldSumToOne(FreqParamGroup *) {
			throw MTHException("Model::FreqParamChangesShouldSumToOne(FreqParamGroup *p) has been called");
		}
	
		virtual unsigned translateLocalStateCodeToGlobal(unsigned localSC) const {
			return localSC;
		}
	protected :
		unsigned addSharedParam(Parameter *p)
		{
			this->allParamVec.push_back(p);
			return allParamVec.size() - 1;
		}
		unsigned addOwnedParam(Parameter *p)
		{
			ownedParamVec.push_back(p);
			return addSharedParam(p);
		}
		unsigned nstates;
		unsigned nparams;
		unsigned nfreeparams;
		unsigned nfreqParamGroups;
		
		std::vector<Parameter *> allParamVec;
		std::vector<Parameter *> ownedParamVec;
		double ***pMatrix; //pointer to the matrices, there will ngamma categories of matrices, essentially a vector, but not using STL to allow pointer math
		bool pmatcalc;
		FreqParamGroup *stateFreqs;
};


//Eigen values stuff pinched from John's Mr. Bayes


class ModelWEig : public Model {
	public :

		ModelWEig(unsigned n);
		ModelWEig(unsigned n, double *EIEMPreAlloc); //memory saving by sharing EIEM
		ModelWEig(unsigned n, double **vecPreAlloc, double ***matPreAlloc, complex ***complexPreAlloc, double ***pMat); //for shared memory	
		ModelWEig(unsigned n, double pi); //some invariant sites
		ModelWEig(unsigned n, unsigned ncats, double gamsh); //gamma rates
		ModelWEig(unsigned n, double pi, unsigned ncats, double gamsh); //gamma +pinvariant
		virtual ~ModelWEig();
		
		virtual void UpdatePMatrix(double **, double);
		virtual void UpdatePMatrix(double **pmats, double blen, unsigned onlycolumn);
		virtual void UpdatePRow(double **pmats, double blen, unsigned onlycolumn);
		
		virtual void CalculateQ()=0;
		virtual void InitializeParameters()=0;
		void SharedConstruction(unsigned n);

		void ParameterHasChanged(Parameter *p) {
			if (p == this->gammashape)
				this->CalculateRates();
			else if (p != this->pinv) {
				this->eigencalc = false;
				this->qmatcalc=false;
			}
		}
		void ParameterHasChanged(FreqParamGroup *) {
			eigencalc = false;
			qmatcalc = false;
		}
	protected:

		complex **CEigenVector, **CInvEigenVector;
		double *REigenValues;
		double *ImEigenValues;
		double **REigenVector;
		double **InvEigenVector;
		double **qMatrix;
		double *EigInvEigMat; //n x n x n matrix of multiplications to cut down on recalcs
		double *EigExp; //n sized vector of exp(Eigenvalues *blen ) to cut down on reallocation of memory everytime blen changes
		bool eigencalc, qmatcalc, IsComplex;
};

} // namespace bull

#endif
