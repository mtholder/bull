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
 

#ifndef SSRFSUBMODCODON
#define SSRFSUBMODCODON

#include <ostream>
#include <fstream>


#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 

#define MINIMUMAAFREQABOVEZERO 1.0e-12

#include "model.hpp"
#include "ssrf_codon_sub_mod.hpp"

namespace bull {

class SSRFCodonSubModSet;
class SSRFCodonSubMod : public ModelWEig
{
	public:
		static bool mutParamsDirty;
		static unsigned geneticCode;
	
		static double defFreqA,defFreqC,defFreqG,defrAC,defrAG,defrAT,defrCG,defrCT,defrGT,
			defA,defC,defD,defE,defF,defG,defH,defI,defK,defL,defM,defN,defP,defQ,defR,defS,defT,
			defV,defW,defY,defStop,defRate;
	
		SSRFCodonSubMod(unsigned ncods, unsigned naas, const bool *codP, const bool *aaP, SSRFCodonSubModSet & ssrfSet, const std::vector<double> & ssrfParams, unsigned nSites);
		virtual ~SSRFCodonSubMod();
		void CalculateQ();
			//written to allow odd coding of SSRFCodonSubModel which overrides
		void EncodeACharacter(short *dest,short *inp,unsigned datatype,bool keepGap) ;
		void FreqParamChangesShouldSumToOne(FreqParamGroup *p);
		unsigned GetEncodingType() const;
		unsigned GetNStates() const;
		virtual const double * const * GetStateFreqsConst() const;
		void InitializeParameters();
		unsigned NumShortsPerCharacter();
		unsigned NumStatesInLastShort();
		void ParameterHasChanged(FreqParamGroup *p) ;
		void ParameterHasChanged(Parameter *p);
		void RecalculateMutMatrix();
		static void set_default(double *allParams);
		inline void UpdatePmat(double b);
		inline void UpdatePmat(double b,unsigned onlycol);
		inline void UpdatePMatrix(double **pm,double b); //for externally supplied Pmatrix
		inline void UpdatePMatrix(double **pm,double b,unsigned onlycol);
		inline void UpdatePRow(double **pm,double b,unsigned onlycol);
		
		
		//AlertSharedMemory() is called because the EigInvEigMult matrix is shared, if the eigen vectors don't need to be recalculated
		//(the mod parameter haven't changed) this will need to be called whenever the model being scored is changed called by Tree::PrepareLikeAttrAndAllNodes(unsigned partnum,unsigned modn)
		void AlertSharedMemory() const;
		double GetMultiplier() const;
		void SetMultiplier(double x);
		void AddAminoAcidFreqs(std::ofstream &dest);
		
		void SetBlenMultBasedOnSBLI(void);
		void Summarize(std::ofstream &outpf);
		void Fill4TaxonSpectrum(double *spec, double bOne, double bTwo, double bInt, double bThr, double bFou);
		void Fill4TaxonSpectrum(double *specOne, double *specTwo, double *specThr, double bOne, double bTwo, double bInt, double bThr, double bFou); //for testing site specific rates
	
		double GetDeviationFromEquilibriumFreq(double branchLen);
		double GetProportionOnSubOptimalHill();
		unsigned GetIndexOfMutationalNeighbor(unsigned inIndex,unsigned numOfMutNeighbor);

		const unsigned nCodonsInProtein; //total # of codons, in the models that share this mutational set of parameters

#		ifdef ALLOWMULTIHITS
			enum PARAM_ENUM {blenMult=0 ,pMultHit,freqA , freqC , freqG , freqT , rAC , rAG , rAT , rCG , rCT , rGT ,fAAs};
#		else
			enum PARAM_ENUM {blenMult=0 ,freqA , freqC , freqG , freqT , rAC , rAG , rAT , rCG , rCT , rGT ,fAAs};
#		endif
	
		virtual unsigned translateLocalStateCodeToGlobal(unsigned localSC) const {
			return codonLocToGlob.at(localSC);
		}
	

	private:

		double CalculateBeta(void);
		void CalculateCodonFreqs();

		const double * baseFreqsAlias[4];
		double *nucRateMatAlias[4][4];
		FreqParamGroup *aaFreqs;
		FreqParamGroup *baseFreqsGroupAlias;
		std::vector<const double *> codFreqPtrs;
		double codonFreqs[64];
		double SSRkContrib,overflowMultiplier;
		Parameter *sharedBrLenInterpreter;
	
		//global codon code: AAA=0 ... TTT=63		global amino acid code A=0 ... Y=19 *=20
		std::vector<unsigned> codonLocToGlob; //codonTrans[i] is the global codon code where i is the local codon code
		std::vector<unsigned> aminoAcidLocToGlob; //aminoAcidTLocToGlob[i] is the global amino acid code where i is the local amino acid code
		unsigned codonGlobToLoc[64]; //thisCodeByStd[i] is the local codon code where i is the global code
		std::vector<unsigned> codLocToAminoAcidLoc; //translationToLocal[i] is the local amino acid code where i is the local codon code

		// in the cpp file there is a global unsigned MitotcodNum[i] that is the global amino acid code when i is the global codon code
		unsigned nPossAAs;
		double **mutMatAlias;
		
		friend class Tree; /*temporary for checkting codLocToAminoAcidLoc to see if there was a change in simulate data function */
};


/*******************************************************************************
* Copies are SHALLOW and destructor does not free.
*
* You must call freeMemory() on only one of the instances that is aliases!
* you can call surrenderThenClear() get rid of all aliases without freeing memory
*/
class SSRFCodonSubModSet
	{
	public:
		SSRFCodonSubModSet()
			:sharedPreAllocated(NULL),
			baseFreqs(NULL),
			mutMatrix(NULL) {
			setMatsToNull();
		}
		
		SSRFCodonSubModSet(const DblVector & gtrParamValues, /* at least 8 doubles long */
						   const DblMatrix & setofSSRFparams, /* full size=nm*21 */
						   const DblVector & multipliers,
						   const double treesc,
						   const unsigned gcode)
			:sharedPreAllocated(NULL),
			baseFreqs(NULL),
			mutMatrix(NULL) {
			setMatsToNull();
			initialize(gtrParamValues, setofSSRFparams, multipliers, treesc, gcode);
		}
		void initialize(const DblVector & gtrParamValues, /* at least 8 doubles long */
						   const DblMatrix & setofSSRFparams, /* full size=nsites (each row has 21 elements) */
						   const DblVector & multipliers,
						   const double treesc,
						   const unsigned gcode);
		
		unsigned getNumAASites() const {
			return modelPtrs.size();
		}
		void freeMemory();

		void surrenderThenClear() {
			modelPtrs.clear();
			sharedParams.clear();
			setMatsToNull();
			sharedPreAllocated = NULL;
			baseFreqs = NULL;
			mutMatrix = NULL;
		}


		double * sharedPreAllocated; // owned (but all models alias this)
	private:
		void setMatsToNull() {
			nucRateMat[0][0] = nucRateMat[0][1] = nucRateMat[0][2] = nucRateMat[0][3] = NULL;
			nucRateMat[1][0] = nucRateMat[1][1] = nucRateMat[1][2] = nucRateMat[1][3] = NULL;
			nucRateMat[2][0] = nucRateMat[2][1] = nucRateMat[2][2] = nucRateMat[2][3] = NULL;
			nucRateMat[3][0] = nucRateMat[3][1] = nucRateMat[3][2] = nucRateMat[3][3] = NULL;
			baseFreqMat[0] = baseFreqMat[1] = baseFreqMat[2] = baseFreqMat[3] = NULL;
		}

		std::vector<SSRFCodonSubMod *> modelPtrs;
		double * nucRateMat[4][4]; //alias
		double * baseFreqMat[4]; // alias
		std::vector<Parameter *> sharedParams; //owned by this set -- shared by the SSRFCodonSubMod "submodels"
		FreqParamGroup * baseFreqs;  //owned
		double ** mutMatrix; // owned psdmatrix

		friend class SSRFCodonSubMod;
		friend class BullKernel;
	};

void ModifyBlenMultSoBranchesAreScaled(SSRFCodonSubMod **modArr,unsigned ns);
unsigned GetNumberOfPossibleAminoAcids(bool *codOfObsAA,bool *possAAs);
unsigned GetNumberOfPossibleCodons(bool *possAAs);
void SetCodeRelatedGlobals(unsigned gcode);
void GetNeutralAminoAcidFreq(double *aafreq,double *baseFreqParams);
//bool CheckNotZero(double x);
unsigned GetGlobalIndexOfMutationalNeighbor(unsigned inIndex,unsigned numOfMutNeighbor);

} //namespace bull {


#endif
