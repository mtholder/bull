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
 

#include <fstream>
#include <vector>

#include "ssrf_codon_sub_mod.hpp"
#include "basic_bull.hpp"
#include "lin_alg.hpp"
#include "matrices.hpp"

using std::cout;
using namespace bull;

double SSRFCodonSubMod::defFreqA(.25);
double SSRFCodonSubMod::defFreqC(.25);
double SSRFCodonSubMod::defFreqG(.25);
double SSRFCodonSubMod::defrAC(1);
double SSRFCodonSubMod::defrAG(2);
double SSRFCodonSubMod::defrAT(1);
double SSRFCodonSubMod::defrCG(1);
double SSRFCodonSubMod::defrCT(2);
double SSRFCodonSubMod::defrGT(1);
double SSRFCodonSubMod::defA(.05);
double SSRFCodonSubMod::defC(.05);
double SSRFCodonSubMod::defD(.05);
double SSRFCodonSubMod::defE(.05);
double SSRFCodonSubMod::defF(.05);
double SSRFCodonSubMod::defG(.05);
double SSRFCodonSubMod::defH(.05);
double SSRFCodonSubMod::defI(.05);
double SSRFCodonSubMod::defK(.05);
double SSRFCodonSubMod::defL(.05);
double SSRFCodonSubMod::defM(.05);
double SSRFCodonSubMod::defN(.05);
double SSRFCodonSubMod::defP(.05);
double SSRFCodonSubMod::defQ(.05);
double SSRFCodonSubMod::defR(.05);
double SSRFCodonSubMod::defS(.05);
double SSRFCodonSubMod::defT(.05);
double SSRFCodonSubMod::defV(.05);
double SSRFCodonSubMod::defW(.05);
double SSRFCodonSubMod::defY(.05);
double SSRFCodonSubMod::defStop(.00);
double SSRFCodonSubMod::defRate(1);
bool   SSRFCodonSubMod::mutParamsDirty(true);
unsigned SSRFCodonSubMod::geneticCode((unsigned) GenCode(MITO)); // TAH changed Mito to MITO

unsigned CodByaa[21][6];
unsigned NCodByAA[21];
unsigned gCodonNumber[64];
unsigned tfirbase[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
unsigned tsecbase[]={0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
unsigned tthibase[]={0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};


void SSRFCodonSubModSet::freeMemory()
{
	for (std::vector<SSRFCodonSubMod *>::iterator mIt = this->modelPtrs.begin(); mIt != this->modelPtrs.end(); ++mIt) {
		SSRFCodonSubMod * p = *mIt;
		delete p;
	}
	this->modelPtrs.clear();
	
	delete [] this->sharedPreAllocated;
	this->sharedPreAllocated = NULL;
	
	delete this->baseFreqs;
	this->baseFreqs = NULL;
	
	this->setMatsToNull();
	
	for (std::vector<Parameter *>::iterator pIt = this->sharedParams.begin(); pIt != this->sharedParams.end(); ++pIt) {
		Parameter * parm = *pIt;
		delete parm;
	}
	sharedParams.clear();
	
	if (this->mutMatrix)
		free_psdmatrix(this->mutMatrix);
	this->mutMatrix = NULL;
}

void SSRFCodonSubModSet::initialize(
	const DblVector & gtrp, /* at least N_MUT_PARAMS doubles long freq(A,C,G) r(AC, AG, AT, CG, CT, GT), pMultHit */
	const DblMatrix & setofSSRFparams, /* full size=nm*21 */
	const DblVector & multipliers,
	const double treesc,
	const unsigned gcode)
{
	const unsigned nAASites = setofSSRFparams.size();
	freeMemory();
	try {
		this->modelPtrs.resize(nAASites, 0L);
		this->sharedPreAllocated = new double [262144]; //64^3
		
		this->sharedParams.resize(SSRFCodonSubMod::rGT + 2, 0L);
		this->sharedParams[SSRFCodonSubMod::blenMult] = new PositiveParameter(treesc,par(MIN)|par(CUR),1.0);
	
#		ifdef ALLOWMULTIHITS
			this->sharedParams[SSRFCodonSubMod::pMultHit] = new PositiveParameter(gtrp[MULTI_HIT_PARAM_INDEX],par(MIN)|par(CUR),0.0);
#		endif

		this->sharedParams[SSRFCodonSubMod::freqA] = new FullParameter(gtrp[0], bull::BULL_SMALL_DBL,1.0-bull::BULL_SMALL_DBL,par(MIN)|par(MAX)|par(CUR),SSRFCodonSubMod::defFreqA,"freqA");
		this->sharedParams[SSRFCodonSubMod::freqC] = new FullParameter(gtrp[1], bull::BULL_SMALL_DBL,1.0-bull::BULL_SMALL_DBL,par(MIN)|par(MAX)|par(CUR),SSRFCodonSubMod::defFreqC,"freqC");
		this->sharedParams[SSRFCodonSubMod::freqG] = new FullParameter(gtrp[2], bull::BULL_SMALL_DBL,1.0-bull::BULL_SMALL_DBL,par(MIN)|par(MAX)|par(CUR),SSRFCodonSubMod::defFreqG,"freqG");
		this->sharedParams[SSRFCodonSubMod::freqT] = new FullParameter(1.0-gtrp[0]-gtrp[1]-gtrp[2],bull::BULL_SMALL_DBL,1.0-bull::BULL_SMALL_DBL,par(MIN)|par(MAX)|par(CUR),1.0-SSRFCodonSubMod::defFreqA-SSRFCodonSubMod::defFreqC-SSRFCodonSubMod::defFreqG,"freqT");
		this->sharedParams[SSRFCodonSubMod::rAC] = new PositiveParameter(gtrp[3],par(MIN)|par(CUR),SSRFCodonSubMod::defrAC);
		this->sharedParams[SSRFCodonSubMod::rAG] = new PositiveParameter(gtrp[4],par(MIN)|par(CUR),SSRFCodonSubMod::defrAG);
		this->sharedParams[SSRFCodonSubMod::rAT] = new PositiveParameter(gtrp[5],par(MIN)|par(CUR),SSRFCodonSubMod::defrAT);
		this->sharedParams[SSRFCodonSubMod::rCG] = new PositiveParameter(gtrp[6],par(MIN)|par(CUR),SSRFCodonSubMod::defrCG);
		this->sharedParams[SSRFCodonSubMod::rCT] = new PositiveParameter(gtrp[7],par(MIN)|par(CUR),SSRFCodonSubMod::defrCT);
		this->sharedParams[SSRFCodonSubMod::rGT] = new PositiveParameter(gtrp[8],par(MIN)|par(CUR),SSRFCodonSubMod::defrGT);
		this->sharedParams[SSRFCodonSubMod::rGT+1] = new PositiveParameter(0.0,par(MIN)|par(CUR),1.0); //the branch length interpreter
	
		this->baseFreqs = new FreqParamGroup(4, &(this->sharedParams[SSRFCodonSubMod::freqA]));
		
		SSRFCodonSubMod::geneticCode = gcode;
		SetCodeRelatedGlobals(gcode);
		SSRFCodonSubMod::mutParamsDirty=true;
		
		this->nucRateMat[0][0] = this->nucRateMat[1][1] = this->nucRateMat[2][2] = this->nucRateMat[3][3] = NULL;
		this->nucRateMat[0][1] = this->nucRateMat[1][0] = &(this->sharedParams[SSRFCodonSubMod::rAC]->val);
		this->nucRateMat[0][2] = this->nucRateMat[2][0] = &(this->sharedParams[SSRFCodonSubMod::rAG]->val);
		this->nucRateMat[0][3] = this->nucRateMat[3][0] = &(this->sharedParams[SSRFCodonSubMod::rAT]->val);
		this->nucRateMat[1][2] = this->nucRateMat[2][1] = &(this->sharedParams[SSRFCodonSubMod::rCG]->val);
		this->nucRateMat[1][3] = this->nucRateMat[3][1] = &(this->sharedParams[SSRFCodonSubMod::rCT]->val);
		this->nucRateMat[2][3] = this->nucRateMat[3][2] = &(this->sharedParams[SSRFCodonSubMod::rGT]->val);
		
		/*baseFreqsAlias = new double[4]; */
		this->baseFreqMat[0] = &(this->sharedParams[SSRFCodonSubMod::freqA]->val);
		this->baseFreqMat[1] = &(this->sharedParams[SSRFCodonSubMod::freqC]->val);
		this->baseFreqMat[2] = &(this->sharedParams[SSRFCodonSubMod::freqG]->val);
		this->baseFreqMat[3] = &(this->sharedParams[SSRFCodonSubMod::freqT]->val);
	
		this->mutMatrix = psdmatrix(64);
		
		bool codOfObsAA[64];
		bool possAAs[21];
		
		for (unsigned i=0; i < 64; i++)
			codOfObsAA[i] = false;
		for (unsigned i=0; i < 21; i++)
			possAAs[i] = false;
	
		for (unsigned i = 0 ; i < nAASites; i++) {
			double tot = 0.0;
			DblVector ssrfParamRow = setofSSRFparams[i];
			for (unsigned j = 0; j < 20; j++)
				tot += ssrfParamRow[j];
	
			for (unsigned k = 0; k < 64; k++) {
				const unsigned nc = gCodonNumber[k];
				if (nc < 20) {
#					ifdef TESTINGINFERENCE
						const bool bv = true; //no packing at all if comparing two different models for the same site (in case on is zero for a state that the other has)
#					else
						const bool bv = (ssrfParamRow[nc] > 0.0);
#					endif
					codOfObsAA[k] = bv;
					possAAs[nc] = bv;
				}
			}
			const bool sbv = (fabs(1.0-tot) >= bull::BULL_SMALL_DBL);
			possAAs[20] = sbv;
			const unsigned * stopCodons = CodByaa[20];
			codOfObsAA[stopCodons[0]] = codOfObsAA[stopCodons[1]] = codOfObsAA[stopCodons[2]] = codOfObsAA[stopCodons[3]] = sbv;
	
			const unsigned nPossibleAAs = GetNumberOfPossibleAminoAcids(codOfObsAA, possAAs);
			const unsigned nPossibleCodons = GetNumberOfPossibleCodons(possAAs);
			this->modelPtrs[i] = new SSRFCodonSubMod(nPossibleCodons, nPossibleAAs, codOfObsAA, possAAs, *this, ssrfParamRow, nAASites);
			if (!multipliers.empty())	
				this->modelPtrs[i]->SetMultiplier(multipliers[i]);
		}
	}
	catch (...)
		{
		this->freeMemory();
		throw;
		}
}



SSRFCodonSubMod::SSRFCodonSubMod(
	unsigned ncods,
	unsigned naas,
	const bool * codP,
	const bool * aaP,
	SSRFCodonSubModSet & ssrfSet,
	const std::vector<double> & ssrfParams,
	unsigned nSites)
	: ModelWEig(ncods, ssrfSet.sharedPreAllocated),
	nCodonsInProtein(nSites),
	aaFreqs(NULL),
	codFreqPtrs(64, (const double *)(NULL))
{
	Parameter ** gtrParams = &(ssrfSet.sharedParams[0]);
	mutMatAlias = ssrfSet.mutMatrix;
	baseFreqsGroupAlias = ssrfSet.baseFreqs;
	const double * const * sf = ssrfSet.baseFreqs->GetStateFreqsConst();
	for (unsigned i = 0; i < 4; ++i)
		baseFreqsAlias[i] = sf[i];
	for (unsigned i = 0; i < 4; ++i)
		for (unsigned j = 0; j < 4; ++j)
			nucRateMatAlias[i][j] = ssrfSet.nucRateMat[i][j];
	for (unsigned i = 0; i < 64; ++i)
		codFreqPtrs[i] = &(codonFreqs[i]);
	assert(sizeof(short) == 2); //if this fails SSRFCodonSubMod::NumShortsPerCharacter will have to be changed
#	ifdef ALLOWMULTIHITS
		nparams = 11 + naas;
#	else
		nparams = 10 + naas;
#	endif
	nPossAAs = naas;
	
	/* The order of parameters must be maintained here!!
	* Must agree with PARAM_ENUM
	*/
	addSharedParam(gtrParams[SSRFCodonSubMod::blenMult]);
#	ifdef ALLOWMULTIHITS
		addSharedParam(gtrParams[SSRFCodonSubMod::pMultHit]);
#	endif
	addSharedParam(gtrParams[SSRFCodonSubMod::freqA]);
	addSharedParam(gtrParams[SSRFCodonSubMod::freqC]);
	addSharedParam(gtrParams[SSRFCodonSubMod::freqG]);
	addSharedParam(gtrParams[SSRFCodonSubMod::freqT]);
	addSharedParam(gtrParams[SSRFCodonSubMod::rAC]);
	addSharedParam(gtrParams[SSRFCodonSubMod::rAG]);
	addSharedParam(gtrParams[SSRFCodonSubMod::rAT]);
	addSharedParam(gtrParams[SSRFCodonSubMod::rCG]);
	addSharedParam(gtrParams[SSRFCodonSubMod::rCT]);
	addSharedParam(gtrParams[SSRFCodonSubMod::rGT]);

	/*
	* End of "The order of parameters must be maintained here!!" section
	**/
	sharedBrLenInterpreter = gtrParams[SSRFCodonSubMod::rGT+1];
	overflowMultiplier = 1.0;
	SSRkContrib = 0.0;
	
	nfreqParamGroups = 2;
	
	codonLocToGlob.assign(nstates, 0);
	codLocToAminoAcidLoc.assign(nstates, 0);
	
	unsigned thisCodonInd=0;
	for (unsigned i = 0; i < 64; i++) {
		if (codP[i]) {
			assert(thisCodonInd < nstates);
			codonGlobToLoc[i] = thisCodonInd;
			codonLocToGlob[thisCodonInd] = i;
			thisCodonInd++;
		}
		else
			codonGlobToLoc[i] = UINT_MAX;
	}
	
	aminoAcidLocToGlob.clear();
	for (unsigned i = 0; i < 21; i++) {
		if (aaP[i])
			aminoAcidLocToGlob.push_back(i);
	}

	for (unsigned i = 0; i < nstates; i++) {
		const unsigned translated = gCodonNumber[codonLocToGlob[i]];
		for (unsigned j = 0; j < nPossAAs; j++) {
			if (translated == aminoAcidLocToGlob[j]) {
				codLocToAminoAcidLoc[i]=j;
				break;
			}
		}
	}

	double OneMinusOthers=1.0;
	for (unsigned j = 0; j < ssrfParams.size(); j++)
		OneMinusOthers -= ssrfParams[j];
	if (nPossAAs == 1)
		addOwnedParam(new FullParameter(1.0, 0.0, 1.0, par(MIN)|par(MAX)|par(CUR)|par(FIX), defA, "Amino"));
	else {
		for (unsigned i=0; i < nPossAAs; i++) {
			FullParameter * nfp;
			if (aminoAcidLocToGlob[i] == 20)
				nfp = new FullParameter(OneMinusOthers,0.0,1.0,par(MIN)|par(MAX)|par(CUR),defA,"Stop");
			else 
				nfp = new FullParameter(ssrfParams[aminoAcidLocToGlob[i]],0.0,1.0,par(MIN)|par(MAX)|par(CUR),defA,"Amino");
			addOwnedParam(nfp);
		}
	}
	
	aaFreqs = new FreqParamGroup(nPossAAs, &(allParamVec[SSRFCodonSubMod::fAAs]));
	
	CalculateQ(); //initializes RateConst parameter and codon freqs
}


void SSRFCodonSubMod::CalculateCodonFreqs()
{	
	double numerators[6];
	for (unsigned i = 0; i< 64; ++i)
		codonFreqs[i] = 0.0;
	for (unsigned thisaa=0; thisaa < nPossAAs; thisaa++) {
		double aaf = allParamVec[SSRFCodonSubMod::fAAs+thisaa]->val;
		unsigned ndegcod=NCodByAA[aminoAcidLocToGlob[thisaa]];
		double denominator=0.0;
		for (unsigned thiscodon=0; thiscodon < ndegcod; thiscodon++) {
			unsigned thisglobcod = CodByaa[aminoAcidLocToGlob[thisaa]][thiscodon];
			assert(thisglobcod!=65);
			numerators[thiscodon] = *(baseFreqsAlias[tfirbase[thisglobcod]]);
			numerators[thiscodon] *= *(baseFreqsAlias[tsecbase[thisglobcod]]);
			numerators[thiscodon] *= *(baseFreqsAlias[tthibase[thisglobcod]]);
			denominator += numerators[thiscodon];
		}
		for (unsigned thiscodon=0; thiscodon < ndegcod; thiscodon++) {
			unsigned thisglobcod = CodByaa[aminoAcidLocToGlob[thisaa]][thiscodon];
			assert(thisglobcod!=65);
			codonFreqs[codonGlobToLoc[thisglobcod]] = aaf*numerators[thiscodon]/denominator;
		}
	}
	double tot=0.0;
	for (unsigned i = 0; i< 64; ++i)
		tot += codonFreqs[i];
	if (!(fabs(tot-1.0)<0.00000001))
		throw MTHException("Codon Frequencies don't add up to 1");	
}

SSRFCodonSubMod::~SSRFCodonSubMod()
{
	delete aaFreqs;
	EigInvEigMat=NULL; //shared memory
}

////////////////////////////////////////////////////////////////////////////////
//	Beta is a factor that is multiplied to the branchlength to make the eqns in 
//	Molecular Systematics work.	 This the constraint that Sum qii * freq(i) = -1 this ensures that the branches are in length
//	that are expected numbers of changes  
//	SSRkContrib is the defined k constant, this function maintains the rates
//	For GTR the rate params are scaled down to maintain their same ratio, but make beta=1.0
////////////////////////////////////////////////////////////////////////////////
double	SSRFCodonSubMod::CalculateBeta()
{	double x = 0.0;
	double *qm = *qMatrix;
	//assumes that calculate code freq; and that the rows sum to zero
	for (unsigned row=0; row <nstates; row++) {
		x -= codonFreqs[row]**qm;
		qm += nstates+1;
	}
	return x;
}

#undef OLDWAYOFMULTIHITS
#define PMUT 0.00001
void SSRFCodonSubMod::RecalculateMutMatrix()
{	register unsigned i,j,k,l,m,n;
	
#	ifdef ALLOWMULTIHITS
#		ifdef OLDWAYOFMULTIHITS
			double ssmm[16];
			for (i=0; i < 4; i++)
				ssmm[5*i]=0.0;
			for (i=0; i < 4; i++)
				for (j=0; j < 4; j++)
					if (i!=j)
						ssmm[5*i]+=ssmm[4*i+j]=*nucRateMatAlias[i][j]**baseFreqsAlias[j];
			double temp=0.0;
			double *tempp;
			for (i=0; i < 4; i++)
				temp+=ssmm[5*i]**baseFreqsAlias[i]; //
			temp=PMUT/temp;
			tempp=ssmm;
			for (i=0; i < 16; i++)
				*tempp++*=temp;
			//ssmm is now scaled so that the weighted average of p(mutation) == PMUT (mutations are rare so a discrete-time motivated producton of multi -site rates is reasonable
			for (i=0; i < 4; i++)
				{ssmm[5*i]=1.0-ssmm[5*i];
				if (ssmm[5*i]<.95)
					throw BadParams("After scaling a nucleotide mutation prob is >.05");
				}
			
			tempp=mutMatAlias[0];
			// now treating allParamVec[SSRFCodonSubMod::pMultHit] as a multiplier to modify the rate of double/ triple hits
			//its affect is dependent on PMUT (which is a define so it won't change in a run) if allParamVec[SSRFCodonSubMod::pMultHit]is
			//chance of one mut is PMUT, double is simply (PMUT^2)* allParamVec[SSRFCodonSubMod::pMultHit] , three is (PMUT^3) * allParamVec[SSRFCodonSubMod::pMultHit]^2, 
			//this is a quasi-probability based method done to mimic AH's implementation (if set to 1.0 you get his version)
			
			for (i=0; i < 4; i++)
			for (j=0; j < 4; j++)
			for (k=0; k < 4; k++)
			for (l=0; l < 4; l++)
			for (m=0; m < 4; m++)
			for (n=0; n < 4; n++)
				{*tempp=ssmm[4*i+l]*ssmm[4*j+m]*ssmm[4*k+n];
				if (i!=l)
					{if(j!=m)
						{if(k!=n)
							*tempp*=(allParamVec[SSRFCodonSubMod::pMultHit]->val)*(allParamVec[SSRFCodonSubMod::pMultHit]->val);
						else	
							*tempp*=(allParamVec[SSRFCodonSubMod::pMultHit]->val);
						}
					else
						{if(k!=n)
							*tempp*=(allParamVec[SSRFCodonSubMod::pMultHit]->val);
						}
					}
				else
					{if(j!=m)
						{if(k!=n)
							*tempp*=(allParamVec[SSRFCodonSubMod::pMultHit]->val);
						}
					}
				tempp++;
				}
			//last step to scale mut rates up to reasonable values (weighted mu is 1.0)
			tempp=mutMatAlias[0];
			temp=0.0;
			unsigned row=0;
			for (i=0; i < 4; i++)
			for (j=0; j < 4; j++)
			for (k=0; k < 4; k++)
				{mutMatAlias[row][row]=0.0;
				for (l=0; l < 64; l++)
					if (l!=row)
						mutMatAlias[row][row]-=mutMatAlias[row][l];
				temp-=mutMatAlias[row][row]**baseFreqsAlias[i]**baseFreqsAlias[j]**baseFreqsAlias[k];
				row++;
				}
			
			for (unsigned i=0; i < 64; i++)
				for (unsigned j=0; j < 64; j++)
					mutMatAlias[i][j]/=temp;
#		else
			
			//New version of multihits, doesn't rely on quasi discrete approximation.  all single mutation events 
			//	occur at 1-PrMultiHit-PrMultiHit^2, and all double hits occur at PrMultiHit, and all triple hits occur at PrMultiHit^2
			
			double *tempp;
			tempp=mutMatAlias[0];
			double pmhsq;
			pmhsq=(allParamVec[SSRFCodonSubMod::pMultHit]->val)*(allParamVec[SSRFCodonSubMod::pMultHit]->val);
			double psingle=1.0-(allParamVec[SSRFCodonSubMod::pMultHit]->val)-pmhsq;
		//	first step to scale mut rates up to consistent meaning to the  MultiHit parameter prob not necessary because matrix is scaled again
			double scalednrm[4][4]; 
			double temp=0.0;
			for (unsigned i=0; i < 4; i++)
				for (j=0; j < 4; j++)
					{if(i!=j)
						{scalednrm[i][j]=*nucRateMatAlias[i][j];
						temp+=scalednrm[i][j]**baseFreqsAlias[j]**baseFreqsAlias[i];
						}
					else
						scalednrm[i][j]=1.0;
					}
			for (unsigned i=0; i < 4; i++)
				for (j=0; j < 4; j++)
					{if(i!=j)
						scalednrm[i][j]/=temp;
					}
			for (i=0; i < 4; i++)
			for (j=0; j < 4; j++)
			for (k=0; k < 4; k++)
			for (l=0; l < 4; l++)
			for (m=0; m < 4; m++)
			for (n=0; n < 4; n++)
				{
				if (i!=l)						//first changed
					{if(j!=m)					//second changed
						{if(k!=n)				// all three changed 
							*tempp++=scalednrm[i][l]**baseFreqsAlias[l]*scalednrm[j][m]**baseFreqsAlias[m]*scalednrm[k][n]**baseFreqsAlias[n]*pmhsq;
						else					//first and second changed
							*tempp++=scalednrm[i][l]**baseFreqsAlias[l]*scalednrm[j][m]**baseFreqsAlias[m]*(allParamVec[SSRFCodonSubMod::pMultHit]->val);
						}
					else
						{if(k!=n)				//first and third changed - happens more or less at the at the three mutation rate
							*tempp++=scalednrm[i][l]**baseFreqsAlias[l]*scalednrm[k][n]**baseFreqsAlias[n]*pmhsq;
						else
							*tempp++=scalednrm[i][l]**baseFreqsAlias[l]*psingle;	//just first changed
						}
					}
				else
					{if(j!=m)					//sec changed
						{if(k!=n)				//sec and third changed
							*tempp++=scalednrm[j][m]**baseFreqsAlias[m]*scalednrm[k][n]**baseFreqsAlias[n]*(allParamVec[SSRFCodonSubMod::pMultHit]->val);
						else					//just sec changed
							*tempp++=scalednrm[j][m]**baseFreqsAlias[m]*psingle;
						}
					else
						{if(k!=n)				//just third changed
							*tempp++=scalednrm[k][n]**baseFreqsAlias[n]*psingle;
						else					//none changed - on diagonal set to zero
							*tempp++=0.0;
						}
					}
				}
		//last step to scale mut rates up to reasonable values (weighted mu is 1.0)
			tempp=mutMatAlias[0];
			temp=0.0;
			unsigned row=0;
			for (i=0; i < 4; i++)
			for (j=0; j < 4; j++)
			for (k=0; k < 4; k++)
				{mutMatAlias[row][row]=0.0;
				for (l=0; l < 64; l++)
					if (l!=row)
						mutMatAlias[row][row]-=mutMatAlias[row][l];
				temp-=mutMatAlias[row][row]**baseFreqsAlias[i]**baseFreqsAlias[j]**baseFreqsAlias[k];
				row++;
				}
			for (unsigned i=0; i < 64; i++)
				for (unsigned j=0; j < 64; j++)
					mutMatAlias[i][j]/=temp;
#		endif
#	else
		double *tempp;
		tempp=mutMatAlias[0];
		for (i=0; i < 4; i++) {
			for (j=0; j < 4; j++) {
				for (k=0; k < 4; k++) {
					for (l=0; l < 4; l++) {
						for (m=0; m < 4; m++) {
							for (n=0; n < 4; n++) {
								double x;
								if (i!=l) 
									x = ((j!=m || k!=n) ? 0.0 : *nucRateMatAlias[i][l] * *baseFreqsAlias[l]);
								else {
									if (j!=m)
										x = (k != n ? 0.0 : *nucRateMatAlias[j][m]**baseFreqsAlias[m]);
									else
										x = (k != n ? *nucRateMatAlias[k][n]**baseFreqsAlias[n] : 0.0);
								}
								*tempp++ = x;
							}
						}
					}
				}
			}
		}	//last step to scale mut rates up to reasonable values (weighted mu is 1.0)
		tempp=mutMatAlias[0];
		double temp=0.0;
		unsigned row=0;
		for (i=0; i < 4; i++)
		for (j=0; j < 4; j++)
		for (k=0; k < 4; k++)
			{mutMatAlias[row][row]=0.0;
			for (l=0; l < 64; l++)
				if (l!=row)
					mutMatAlias[row][row]-=mutMatAlias[row][l];
			temp-=mutMatAlias[row][row]**baseFreqsAlias[i]**baseFreqsAlias[j]**baseFreqsAlias[k];
			row++;
			}
		for (unsigned i=0; i < 64; i++)
			for (unsigned j=0; j < 64; j++)
				mutMatAlias[i][j]/=temp;
					
#	endif
	mutParamsDirty=false;
}
	
void SSRFCodonSubMod::CalculateQ(void)
{	
	sharedBrLenInterpreter->val-=SSRkContrib/nCodonsInProtein; //subtract out contribution of old qmatrix to mean rate of change
	
	double *qm,tempx;
	CalculateCodonFreqs();
	if (mutParamsDirty)
		RecalculateMutMatrix();
	
	unsigned nsq=nstates*nstates;
	qm=qMatrix[0];
	for ( unsigned i=0; i < nsq; i++)
		*qm++=0.0;
	unsigned currFromFirBase,currFromSecBase,currFromThiBase;
	unsigned globTocod, globFromcod, locTocod, locToAA, locFromAA;

	unsigned currToBase;
	for (unsigned row=0; row < nstates; row++) {
		qMatrix[row][row]=0.0;
		locFromAA=codLocToAminoAcidLoc[row];
		if (allParamVec[SSRFCodonSubMod::fAAs+locFromAA]->val>0.0) {
			globFromcod=codonLocToGlob[row];
			currFromFirBase=tfirbase[globFromcod];
			currFromSecBase=tsecbase[globFromcod];
			currFromThiBase=tthibase[globFromcod];
			
			for (currToBase=0; currToBase < 4; currToBase++) {
				if (currToBase != currFromFirBase) {
					globTocod = 16*currToBase + 4*currFromSecBase + currFromThiBase;
					locTocod = codonGlobToLoc[globTocod];
					if (locTocod != UINT_MAX) {
						locToAA=codLocToAminoAcidLoc[locTocod];
						if (locToAA == locFromAA)
							qMatrix[row][locTocod] = mutMatAlias[globFromcod][globTocod];
						else if (codonFreqs[locTocod]<bull::BULL_SMALL_DBL)
							qMatrix[row][locTocod] = 0.0;
						else {
							tempx = codonFreqs[row]*mutMatAlias[globFromcod][globTocod]/(codonFreqs[locTocod]*mutMatAlias[globTocod][globFromcod]);
							const double sdmult = (fabs(tempx-1.0) < bull::BULL_SMALL_DBL ? 1.0 : log(1.0/tempx)/(1.0-tempx));
							qMatrix[row][locTocod] = mutMatAlias[globFromcod][globTocod] * sdmult;
						}
						qMatrix[row][row] -= qMatrix[row][locTocod];
					}
				}
			}	
			for (currToBase=0; currToBase < 4; currToBase++) {
				if (currToBase != currFromSecBase) {
					globTocod = 16*currFromFirBase + 4*currToBase + currFromThiBase;
					locTocod=codonGlobToLoc[globTocod];
					if (locTocod!=UINT_MAX)	{
						locToAA=codLocToAminoAcidLoc[locTocod];
						if (locToAA == locFromAA) 
							qMatrix[row][locTocod]=mutMatAlias[globFromcod][globTocod];
						else if (codonFreqs[locTocod]<bull::BULL_SMALL_DBL)
							qMatrix[row][locTocod]=0.0;
						else {
							tempx=codonFreqs[row]*mutMatAlias[globFromcod][globTocod]/(codonFreqs[locTocod]*mutMatAlias[globTocod][globFromcod]);
							const double sdmult = (fabs(tempx-1.0) < bull::BULL_SMALL_DBL ? 1.0 : log(1.0/tempx)/(1.0-tempx));
							qMatrix[row][locTocod] = mutMatAlias[globFromcod][globTocod] * sdmult;
						}
						qMatrix[row][row] -= qMatrix[row][locTocod];
					}
				}
			}
			for (currToBase=0; currToBase < 4; currToBase++) {
				if (currToBase!=currFromThiBase) {
					globTocod= 16*currFromFirBase + 4*currFromSecBase + currToBase;
					locTocod=codonGlobToLoc[globTocod];
					if (locTocod != UINT_MAX)	{
						locToAA=codLocToAminoAcidLoc[locTocod];
						if (locToAA == locFromAA)
							qMatrix[row][locTocod] = mutMatAlias[globFromcod][globTocod];
						else if (codonFreqs[locTocod]<bull::BULL_SMALL_DBL)
							qMatrix[row][locTocod] = 0.0;
						else {
							tempx = codonFreqs[row]*mutMatAlias[globFromcod][globTocod]/(codonFreqs[locTocod]*mutMatAlias[globTocod][globFromcod]);
							const double sdmult = (fabs(tempx-1.0) < bull::BULL_SMALL_DBL ? 1.0 : log(1.0/tempx)/(1.0-tempx));
							qMatrix[row][locTocod] = mutMatAlias[globFromcod][globTocod]*sdmult;
						}
						qMatrix[row][row] -= qMatrix[row][locTocod];
					}
				}
			}
		}
		else {
			for (unsigned col=0; col < nstates; col++)
				qMatrix[row][col]=0.0;
		}
	}
	SSRkContrib=CalculateBeta();
	sharedBrLenInterpreter->val+=SSRkContrib/nCodonsInProtein; //add in the contribution of the new qmatrix to mean rate of change
}

void SSRFCodonSubMod::InitializeParameters()
{	RateManager::InitializeParameters();
	for (unsigned i=0; i < nparams; i++)
		{if(allParamVec[i]->StartWithCurrent()) ;
		else if (allParamVec[i]->StartWithRandom())
					throw IncompleteModel("Random Function to initialize parameters isn't available yet");
				//allParamVec[i]->val=SomeRandomNumberFunction();
		else if (allParamVec[i]->StartWithApproximation())
					throw IncompleteModel("Initial approximation of parameters isn't available yet");
		else if (allParamVec[i]->StartWithDefault())
					allParamVec[i]->SetToDefault();
		else	throw BadSettings("No starting value of a parameter has been defined");
		}
}
void SSRFCodonSubMod::EncodeACharacter(short *dest,short *raw,unsigned datatype,bool keepGap) 
{	assert((datatype == EncodingType(MitoCodons) || datatype == EncodingType(NucCodons)) && !keepGap);
		//assumes that a codon code will take up 4 shorts;
	assert(sizeof(short) == 2);
	
	short *temp;
	temp=dest;
	
	*dest=0;
	short smask,fb,sb,tb;
	for (unsigned i=0; i < nstates; i++)
		{fb=sb=tb=1;
		fb <<= tfirbase[codonLocToGlob[i]];
		sb <<= tsecbase[codonLocToGlob[i]];
		tb <<= tthibase[codonLocToGlob[i]];
		if (i%16 == 0)
			{smask=1;
			if (i)
				*++dest=0;
			}
		else	smask <<= 1;
		if ((fb&raw[0]) && (sb&raw[1]) && (tb&raw[2])) 
			(*dest)|=smask;
		}
	bool atLeastOne=false;
	for (unsigned i=0; (!atLeastOne && i < nstates); i++)
		if (temp[i])
			atLeastOne=true;
	if (!atLeastOne)
		temp[0]=0;
}

void SSRFCodonSubMod::UpdatePmat(double b)	
{	ModelWEig::UpdatePMatrix(*pMatrix,allParamVec[SSRFCodonSubMod::blenMult]->val*b);
}

void SSRFCodonSubMod::UpdatePmat(double b,unsigned onlycol)	
{	
	Parameter * blen_multiplier = allParamVec[SSRFCodonSubMod::blenMult];
	assert(blen_multiplier);
	ModelWEig::UpdatePMatrix(*pMatrix,blen_multiplier->val*b,onlycol);
}

void SSRFCodonSubMod::UpdatePMatrix(double **pm,double b)	
{	ModelWEig::UpdatePMatrix(pm,allParamVec[SSRFCodonSubMod::blenMult]->val*b);
}

void SSRFCodonSubMod::UpdatePMatrix(double **pm,double b,unsigned onlycol)	
{	ModelWEig::UpdatePMatrix(pm,allParamVec[SSRFCodonSubMod::blenMult]->val*b,onlycol);
}

void SSRFCodonSubMod::UpdatePRow(double **pm,double b,unsigned onlycol)	
{	ModelWEig::UpdatePRow(pm,allParamVec[SSRFCodonSubMod::blenMult]->val*b,onlycol);
}

unsigned SSRFCodonSubMod::GetEncodingType() const {	
	if (geneticCode == GenCode(MITO)) // TAH changed Mito to MITO
		return EncodingType(MitoCodons);
	assert(geneticCode == GenCode(NUCLEAR)); // TAH changed Nuclear to NUCLEAR
	return EncodingType(NucCodons);
}

unsigned SSRFCodonSubMod::GetNStates() const
{	return nstates;
}

void SSRFCodonSubMod::AlertSharedMemory() const
{	
	if (eigencalc)	{
#		ifdef	CONDENSE_MATRICES
			CalculateAndCondenseEigInvEigMult(*REigenVector,InvEigenVector,nstates,EigInvEigMat);
#		else
#			error "linalg assumes that CONDENSE_MATRICES is 1
			CalculateGlobalEigInvEigMult(*REigenVector,InvEigenVector,nstates,EigInvEigMat);
#		endif
	}
}

double SSRFCodonSubMod::GetMultiplier() const {
	return overflowMultiplier;
}

void SSRFCodonSubMod::SetMultiplier(double x)	
{	overflowMultiplier=x;
}

void SSRFCodonSubMod::ParameterHasChanged(FreqParamGroup *p)	
{
	eigencalc = qmatcalc = false;
	if (baseFreqsGroupAlias == p)
		mutParamsDirty = true;
}

void SSRFCodonSubMod::FreqParamChangesShouldSumToOne(FreqParamGroup *p) 
{
	eigencalc=qmatcalc=false;

	if (baseFreqsGroupAlias == p)
		mutParamsDirty = true;
}

void SSRFCodonSubMod::ParameterHasChanged(Parameter *p) 
{if(p == gammashape)	CalculateRates();
else 
	{if(p!=pinv && p!=allParamVec[SSRFCodonSubMod::blenMult]) 
		eigencalc=qmatcalc=false;
#ifdef ALLOWMULTIHITS
	for (unsigned i=SSRFCodonSubMod::pMultHit; ((!mutParamsDirty) && i < SSRFCodonSubMod::fAAs); i++)
#else
	for (unsigned i=SSRFCodonSubMod::freqA; ((!mutParamsDirty) && i < SSRFCodonSubMod::fAAs); i++)
#endif
		if (p == allParamVec[i])
			mutParamsDirty=true;
	}
}
	
const double * const *SSRFCodonSubMod::GetStateFreqsConst() const
{	
	return &codFreqPtrs[0];
}

unsigned SSRFCodonSubMod::NumShortsPerCharacter()	
{	return (1+(unsigned)floor(((double)nstates-.1)/16.0));
}//written to allow odd coding of SSRFCodonSubModel which overrides

unsigned SSRFCodonSubMod::NumStatesInLastShort() 
{	return 1+((nstates-1)%16);
}//written to allow odd coding of SSRFCodonSubModel which overrides

void SSRFCodonSubMod::SetBlenMultBasedOnSBLI(void)
{	
	allParamVec[SSRFCodonSubMod::blenMult]->val=(3.0/sharedBrLenInterpreter->val);
}

void SSRFCodonSubMod::Summarize(std::ofstream &outpf)
{	
	CalculateQ();
	outpf<<nstates<<"\t";
	outpf<<SSRkContrib<<"\t";
	double rateOne,rateTwo,rateThree,temprOn,temprTw,temprTh,synR,nonSynR,tempSynR,tempNonSynR;
	unsigned currFromFirBase,currFromSecBase,currFromThiBase;
	unsigned globTocod,globFromcod,locTocod,locToAA,locFromAA;
	double pctTimeGTR;
	pctTimeGTR=rateOne=rateTwo=rateThree=synR=nonSynR=0.0;
#ifdef ALLOWMULTIHITS
	assert(0); //unwritten
#endif
	unsigned currToBase;
	for (unsigned row=0; row < nstates; row++)
		{locFromAA=codLocToAminoAcidLoc[row];
		temprOn=temprTw=temprTh=tempSynR=tempNonSynR=0.0;
		if (allParamVec[SSRFCodonSubMod::fAAs+locFromAA]->val>0.0)
			{bool fourfolddeg=true;
			globFromcod=codonLocToGlob[row];
			currFromFirBase=tfirbase[globFromcod];
			currFromSecBase=tsecbase[globFromcod];
			currFromThiBase=tthibase[globFromcod];
			for (currToBase=0; currToBase < 4; currToBase++)//try mutating first base
				if (currToBase!=currFromFirBase)		//only check if it is a mutant
					{globTocod=16*currToBase + 4*currFromSecBase + currFromThiBase;
					locTocod=codonGlobToLoc[globTocod];
					if (locTocod!=UINT_MAX)				//if amino acid could have freq >0.0
						{locToAA=codLocToAminoAcidLoc[locTocod];
						temprOn+=qMatrix[row][locTocod];
						if (locToAA == locFromAA)			//same amino acid
							tempSynR+=qMatrix[row][locTocod];
						else
							tempNonSynR+=qMatrix[row][locTocod];
						}
					}
			for (currToBase=0; currToBase < 4; currToBase++)//try mutating second base
				if (currToBase!=currFromSecBase)		//only check if it is a mutant
					{globTocod=16*currFromFirBase + 4*currToBase + currFromThiBase;
					locTocod=codonGlobToLoc[globTocod];
					if (locTocod!=UINT_MAX)				//if amino acid could have freq >0.0
						{locToAA=codLocToAminoAcidLoc[locTocod];
						temprTw+=qMatrix[row][locTocod];
						if (locToAA == locFromAA)			//same amino acid
							tempSynR+=qMatrix[row][locTocod];
						else
							tempNonSynR+=qMatrix[row][locTocod];
						}
					}
					
			
			for (currToBase=0; currToBase < 4; currToBase++)//try mutating third base
				if (currToBase!=currFromThiBase)		//only check if it is a mutant
					{globTocod=16*currFromFirBase + 4*currFromSecBase + currToBase;
					locTocod=codonGlobToLoc[globTocod];
					if (locTocod!=UINT_MAX)				//if amino acid could have freq >0.0
						{locToAA=codLocToAminoAcidLoc[locTocod];
						temprTh+=qMatrix[row][locTocod];
						if (locToAA == locFromAA)			//same amino acid
							tempSynR+=qMatrix[row][locTocod];
						else
							{fourfolddeg=false;
							tempNonSynR+=qMatrix[row][locTocod];
							}
						}
					}
			
			rateOne+=temprOn*(codonFreqs[row]);
			rateTwo+=temprTw*(codonFreqs[row]);
			rateThree+=temprTh*(codonFreqs[row]);
			synR+=tempSynR*(codonFreqs[row]);
			nonSynR+=tempNonSynR*(codonFreqs[row]);
			if (fourfolddeg)
				pctTimeGTR+=(codonFreqs[row]);
			}
		}

	double varrateOne,varrateTwo,varrateThree,varsynR,varnonSynR;
	varrateOne=varrateTwo=varrateThree=varsynR=varnonSynR=0.0;
#ifdef ALLOWMULTIHITS
	assert(0); //unwritten
#endif
	for (unsigned row=0; row < nstates; row++)
		{locFromAA=codLocToAminoAcidLoc[row];
		double vartemprOn,vartemprTw,vartemprTh,vartempSynR,vartempNonSynR;
		vartemprOn=vartemprTw=vartemprTh=vartempSynR=vartempNonSynR=0.0;
		if (allParamVec[SSRFCodonSubMod::fAAs+locFromAA]->val>0.0)
			{bool fourfolddeg=true;
			globFromcod=codonLocToGlob[row];
			currFromFirBase=tfirbase[globFromcod];
			currFromSecBase=tsecbase[globFromcod];
			currFromThiBase=tthibase[globFromcod];
			for (currToBase=0; currToBase < 4; currToBase++)//try mutating first base
				if (currToBase!=currFromFirBase)		//only check if it is a mutant
					{globTocod=16*currToBase + 4*currFromSecBase + currFromThiBase;
					locTocod=codonGlobToLoc[globTocod];
					if (locTocod!=UINT_MAX)				//if amino acid could have freq >0.0
						{locToAA=codLocToAminoAcidLoc[locTocod];
						vartemprOn+=qMatrix[row][locTocod];
						if (locToAA == locFromAA)			//same amino acid
							vartempSynR+=qMatrix[row][locTocod];
						else
							vartempNonSynR+=qMatrix[row][locTocod];
						}
					}
			for (currToBase=0; currToBase < 4; currToBase++)//try mutating second base
				if (currToBase!=currFromSecBase)		//only check if it is a mutant
					{globTocod=16*currFromFirBase + 4*currToBase + currFromThiBase;
					locTocod=codonGlobToLoc[globTocod];
					if (locTocod!=UINT_MAX)				//if amino acid could have freq >0.0
						{locToAA=codLocToAminoAcidLoc[locTocod];
						vartemprTw+=qMatrix[row][locTocod];
						if (locToAA == locFromAA)			//same amino acid
							vartempSynR+=qMatrix[row][locTocod];
						else
							vartempNonSynR+=qMatrix[row][locTocod];
						}
					}
					
			
			for (currToBase=0; currToBase < 4; currToBase++)//try mutating third base
				if (currToBase!=currFromThiBase)		//only check if it is a mutant
					{globTocod=16*currFromFirBase + 4*currFromSecBase + currToBase;
					locTocod=codonGlobToLoc[globTocod];
					if (locTocod!=UINT_MAX)				//if amino acid could have freq >0.0
						{locToAA=codLocToAminoAcidLoc[locTocod];
						vartemprTh+=qMatrix[row][locTocod];
						if (locToAA == locFromAA)			//same amino acid
							vartempSynR+=qMatrix[row][locTocod];
						else
							{fourfolddeg=false;
							vartempNonSynR+=qMatrix[row][locTocod];
							}
						}
					}
			
			varrateOne+=(vartemprOn-rateOne)*(vartemprOn-rateOne)*(codonFreqs[row]);
			varrateTwo+=(vartemprTw-rateTwo)*(vartemprTw-rateTwo)*(codonFreqs[row]);
			varrateThree+=(vartemprTh-rateThree)*(vartemprTh-rateThree)*(codonFreqs[row]);
			varsynR+=(vartempSynR-synR)*(vartempSynR-synR)*(codonFreqs[row]);
			varnonSynR+=(vartempNonSynR-nonSynR)*(vartempNonSynR-nonSynR)*(codonFreqs[row]);
			}
		
		}
		
	outpf<<rateOne<<"\t"<<sqrt(varrateOne)<<"\t";
	outpf<<rateTwo<<"\t"<<sqrt(varrateTwo)<<"\t";
	outpf<<rateThree<<"\t"<<sqrt(varrateThree)<<"\t";
	
	outpf<<synR<<"\t"<<sqrt(varsynR)<<"\t";
	outpf<<nonSynR<<"\t"<<sqrt(varnonSynR)<<"\t";
	outpf<<pctTimeGTR<<"\t";
	outpf<<"\n";
}


double SSRFCodonSubMod::GetDeviationFromEquilibriumFreq(double branchLen)
{	//this function returns the weighted average of deviation each row of the pmatrix from the equilibrium freqs of each codon.
	// the deviation is the sum of the abs val of the difference in the prob of a state (the col) and its eq. freq.	 the weighting is the 
	// eq freq of that codon (the row)
	// the motivation is to describe how long a branch has to be for this site to reflect the true freqs of amino acids so
	// that the relative difficulty of the problem of inferring parameters can be compared from site to site (to assess which 
	// sites' parameter estimate are likely to be most problematic.
	UpdatePmat(branchLen);
	double dev, tempdev;
	dev=0.0;
	for (unsigned i=0; i < nstates; i++)
		{tempdev=0.0;
		for (unsigned j=0; j < nstates; j++)
			tempdev+=fabs(pMatrix[0][i][j] - codonFreqs[j]);
			
		dev+=tempdev*codonFreqs[i];
		}
	return dev;
}

double	SSRFCodonSubMod::GetProportionOnSubOptimalHill()
{	//return the equ freqs of states which are only on sub optimal hill(s)
	// a hill a collection states whose eqfreqs are less than the eq freq of the peak
	// and have a path of single mutational steps of increasing eq freq to the peak
	
	//find highest peak
	CalculateCodonFreqs();
	double highesEF=0.0;
	unsigned highIndex=0;
	for (unsigned i=0; i < nstates; i++)
		if (codonFreqs[i] > highesEF)
			{highesEF= codonFreqs[i];
			highIndex=i;
			}
	bool onOptPeak[64];
	for (unsigned i=0; i < nstates; i++)
		onOptPeak[i]=false;
	onOptPeak[highIndex]=true;
	std::vector < unsigned> indecesOnOptPeak;
	indecesOnOptPeak.push_back(highIndex);
	
	for (unsigned alreadyOnPeak=0; alreadyOnPeak < indecesOnOptPeak.size(); alreadyOnPeak++) {
		//go through all of the mutational neighbors and add them if their eqFreq is lower
		for (unsigned i=0; i < 9; i++) {
			unsigned neighborstate = GetIndexOfMutationalNeighbor(indecesOnOptPeak[alreadyOnPeak],i);
			if (neighborstate != UINT_MAX && codonFreqs[neighborstate] < codonFreqs[indecesOnOptPeak[alreadyOnPeak]]) {
					//neighbor is lower on hill
					indecesOnOptPeak.push_back(neighborstate);
					onOptPeak[neighborstate]=true;
				}
		}
	}
	double subOptPro=0.0;
	for (unsigned i=0; i < nstates; i++)
		if (!onOptPeak[i])
			subOptPro += codonFreqs[i];
	return subOptPro;
}
unsigned SSRFCodonSubMod::GetIndexOfMutationalNeighbor(unsigned inIndex,unsigned numOfMutNeighbor)
{	
	unsigned startingCod,endingCod;
	startingCod = codonLocToGlob[inIndex];
	endingCod = codonGlobToLoc[GetGlobalIndexOfMutationalNeighbor(startingCod,numOfMutNeighbor)];
	return endingCod;
}

void bull::ModifyBlenMultSoBranchesAreScaled(SSRFCodonSubMod **modArr,unsigned ns)
{		//i've been sloppy about calling CalculateQ, so make sure all of the models are updated
	for (unsigned i=0; i < ns; i++)
		modArr[i]->CalculateQ();
		//assumes that all of the models share the same blenMult and sharedBrLenInterpreter
	(*modArr)->SetBlenMultBasedOnSBLI();
		/*
		ofstream modDes("ModelDescription");
		for (unsigned i=0; i < ns; i++)
			modArr[i]->Summarize(modDes);
			
		modDes.close();
		*/
}

	
unsigned bull::GetNumberOfPossibleAminoAcids(bool *codOfObsAA, bool *possAAs) 
{
	unsigned naas=0;
	for (unsigned i=0; i < 21; i++) {
		if (possAAs[i])
			naas++;
	}
	unsigned prevnaas=UINT_MAX;
	while (prevnaas != naas) {
		prevnaas = naas;
		unsigned codn=0;
		for (unsigned fb=0; fb < 4; fb++) {
			for (unsigned sb=0; sb < 4; sb++) {
				for (unsigned tb=0; tb < 4; tb++) {
					if (!codOfObsAA[codn]) {
						bool mfir=false;
						for (unsigned mutatedbase=0; (!mfir && mutatedbase < 4); mutatedbase++)
							if (codOfObsAA[16*mutatedbase+4*sb+tb])
								mfir=true;
						bool msec=false;
						for (unsigned mutatedbase=0; (!msec && mutatedbase < 4); mutatedbase++)
							if (codOfObsAA[16*fb+4*mutatedbase+tb])
								msec=true;
						bool mthi=false;
						for (unsigned mutatedbase=0; (!mthi && mutatedbase < 4); mutatedbase++)
							if (codOfObsAA[16*fb+4*sb+mutatedbase])
								mthi=true;
						if ((mfir&&msec) ||(mthi&&msec) ||(mfir&&mthi)) {
							naas++;
							possAAs[gCodonNumber[codn]]=true;
							for (unsigned k=0; k < 64; k++)
								if (gCodonNumber[codn] == gCodonNumber[k])
									codOfObsAA[k]=true;
						}
					}
					codn++;
				}
			}
		}
	}
	return naas;
}





unsigned bull::GetNumberOfPossibleCodons(bool *possAAs)
{	unsigned nc=0;
	for (unsigned i=0; i < 21; i++)
		if (possAAs[i])
			nc+=NCodByAA[i];
	return nc;
}

void SSRFCodonSubMod::AddAminoAcidFreqs(std::ofstream &dest)
{	//for printing out amino acid freqs in log file
	double outpaaf[20];
	for (unsigned i=0; i < 20; i++)
		outpaaf[i]=0.0;
	for (unsigned i=0; i < nPossAAs; i++)
		if (aminoAcidLocToGlob[i]<20)
			outpaaf[aminoAcidLocToGlob[i]]=allParamVec[SSRFCodonSubMod::fAAs+i]->val;
	for (unsigned i=0; i < 20; i++)
		{dest<<outpaaf[i]<<" ";
		}
}

void bull::SetCodeRelatedGlobals(unsigned gCode)
{	int **CEcodbyaa;
	int *CENCodByAA;
	int *CECodNum;
	if (gCode == (unsigned) GenCode(MITO)) {  
		CEcodbyaa = GetCodByAA(EncodingType(MitoCodons));
		CENCodByAA = GetNCodByAA(EncodingType(MitoCodons));
		CECodNum = GetCodNum(EncodingType(MitoCodons));
	}
	else if (gCode == (unsigned) GenCode(NUCLEAR)) {
		CEcodbyaa=GetCodByAA(EncodingType(NucCodons));
		CENCodByAA=GetNCodByAA(EncodingType(NucCodons));
		CECodNum=GetCodNum(EncodingType(NucCodons));
	}
	
	else assert(0);
	for (unsigned i=0; i < 21; i++)
		for (unsigned j=0; j < 6; j++)
			CodByaa[i][j]=CEcodbyaa[i][j];
	for (unsigned i=0; i < 21; i++)
			NCodByAA[i]=CENCodByAA[i];
	for (unsigned i=0; i < 64; i++)
		gCodonNumber[i]=CECodNum[i];
	
	/*if(codeToUse == GenCode(Mito))
		{unsigned MitocodByaa[21][6]={
			{36,37,38,39,65,65},
			{57,59,65,65,65,65},
			{33,35,65,65,65,65},
			{32,34,65,65,65,65},
			{61,63,65,65,65,65},
			{40,41,42,43,65,65},
			{17,19,65,65,65,65},
			{13,15,65,65,65,65},
			{0 ,2 ,65,65,65,65},
			{28,29,30,31,60,62},
			{12,14,65,65,65,65},
			{1 ,3 ,65,65,65,65},
			{20,21,22,23,65,65},
			{16,18,65,65,65,65},
			{24,25,26,27,65,65},
			{9 ,11,52,53,54,55},
			{4 ,5 ,6 ,7 ,65,65},
			{44,45,46,47,65,65},
			{56,58,65,65,65,65},
			{49,51,65,65,65,65},
			{8 ,10,48,50,65,65}
			};
		unsigned i,j;
		for (i=0; i < 21; i++)
			for (j=0; j < 6; j++)
				CodByaa[i][j]=MitocodByaa[i][j];
				
		unsigned MitoNCodByAA[]={4,2,2,2,2,4,2,2,2,6,2,2,4,2,4,6,4,4,2,2,4};
		for (i=0; i < 21; i++)
			NCodByAA[i]=MitoNCodByAA[i];
		unsigned MitoCodNum[]=
		{8,11,8,11,16,16,16,16,20,15,20,15,10,7,10,7,13,6,13,6,12,12,12,12,
		14,14,14,14,9,9,9,9,3,2,3,2,0,0,0,0,5,5,5,5,17,17,17,17,20,19,20,19,
		15,15,15,15,18,1,18,1,9,4,9,4};
		
		unsigned mitofirbase[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
		unsigned mitosecbase[]={0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
		unsigned mitothibase[]={0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
		for (i=0; i < 64; i++)
			{gCodonNumber[i]=MitoCodNum[i];
			tfirbase[i]=mitofirbase[i];
			tsecbase[i]=mitosecbase[i];
			tthibase[i]=mitothibase[i];
			}
		
		}
	else
		{assert(GenCode(Nuclear) == codeToUse);
		unsigned NuccodByaa[21][6]={
		{36,37,38,39,65,65},
		{57,59,65,65,65,65},
		{33,35,65,65,65,65},
		{32,34,65,65,65,65},
		{61,63,65,65,65,65},
		{40,41,42,43,65,65},
		{17,19,65,65,65,65},
		{12,13,15,65,65,65},
		{0,2,65,65,65,65},
		{28,29,30,31,60,62},
		{14,65,65,65,65,65},
		{1,3,65,65,65,65},
		{20,21,22,23,65,65},
		{16,18,65,65,65,65},
		{8,10,24,25,26,27},
		{9,11,52,53,54,55},
		{4,5,6,7,65,65},
		{44,45,46,47,65,65},
		{58,65,65,65,65,65},
		{49,51,65,65,65,65},
		{48,50,56,65,65,65}
		};

		unsigned i,j;
		for (i=0; i < 21; i++)
			for (j=0; j < 6; j++)
				CodByaa[i][j]= NuccodByaa[i][j];
				
		unsigned NucNCodByAA[]={4,2,2,2,2,4,2,3,2,6,1,2,4,2,6,6,4,4,1,2,3};
		for (i=0; i < 21; i++)
			NCodByAA[i]= NucNCodByAA[i];
		unsigned NucCodNum[]=
		{8,11,8,11,16,16,16,16,14,15,14,15,7,7,10,7,13,6,13,6,12,12,12,12,
		14,14,14,14,9,9,9,9,3,2,3,2,0,0,0,0,5,5,5,5,17,17,17,17,20,19,20,19,
		15,15,15,15,20,1,18,1,9,4,9,4};
		
		
		unsigned Nucfirbase[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
		unsigned Nucsecbase[]={0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3,0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3};
		unsigned Nucthibase[]={0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
		for (i=0; i < 64; i++)
			{gCodonNumber[i]= NucCodNum[i];
			tfirbase[i]= Nucfirbase[i];
			tsecbase[i]= Nucsecbase[i];
			tthibase[i]= Nucthibase[i];
			}
		
		}
	*/
}

void bull::GetNeutralAminoAcidFreq(double *aafreq,double *baseFreqParams)
{//use only first three base freqs!!
	double baseFr[4];
	baseFr[0]=baseFreqParams[0];
	baseFr[1]=baseFreqParams[1];
	baseFr[2]=baseFreqParams[2];
	baseFr[3]=1-baseFr[2]-baseFr[1]-baseFr[0];
	for (unsigned i=0; i < 21; i++)
		{aafreq[i]=baseFr[tfirbase[CodByaa[i][0]]] * baseFr[tsecbase[CodByaa[i][0]]] * baseFr[tthibase[CodByaa[i][0]]] ; 
		for (unsigned j=1; j < NCodByAA[i]; j++)
			aafreq[i]+=baseFr[tfirbase[CodByaa[i][j]]] * baseFr[tsecbase[CodByaa[i][j]]] * baseFr[tthibase[CodByaa[i][j]]] ; 
		}
	
}

unsigned bull::GetGlobalIndexOfMutationalNeighbor(unsigned inIndex,unsigned numOfMutNeighbor)
{
	assert(numOfMutNeighbor < 9);
//takes and returns global codon num
	unsigned mutbase=0;
	int basesMoved=-1;
	if (numOfMutNeighbor < 3) {
		//mutating third
		while (basesMoved < (int)numOfMutNeighbor) {
			if(mutbase != tthibase[inIndex])
				basesMoved++;
			mutbase++;
		}
		mutbase--;
		return (16*tfirbase[inIndex] + 4*tsecbase[inIndex] + mutbase);
	}
	else if (numOfMutNeighbor < 6) {
		//mutating second
		numOfMutNeighbor-=3;
		while (basesMoved < (int)numOfMutNeighbor) {
			if(mutbase!=tsecbase[inIndex])
				basesMoved++;
			mutbase++;
		}
		mutbase--;
		return (16*tfirbase[inIndex] + 4*mutbase + tthibase[inIndex]);
	}
	numOfMutNeighbor-=6;
	while (basesMoved < (int)numOfMutNeighbor) {
		if(mutbase!=tfirbase[inIndex])
			basesMoved++;
		mutbase++;
	}
	mutbase--;
	return 16*mutbase + 4*tsecbase[inIndex] + tthibase[inIndex];
}
