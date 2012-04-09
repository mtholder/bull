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

#include "bull_kernel.hpp"

#include "bull_parsers.hpp"
#include "tree.hpp"
#include "string_extensions.hpp"
#include "xbull.hpp"


#include "ncl/nxsstring.h"

using namespace std;
using namespace bull;

BullKernel::BullKernel() 
  :charsOpen(false),
  modelsOpen(false),
  taxaOpen(false),
  treesOpen(false)
{
	Reset();
	notifyListeners(BullListener::kInitialized);
}


Tree * BullKernel::FindTreeFromName(std::string s)
{
	for (std::list<Tree *>::iterator trIt = treeList.begin(); trIt != treeList.end(); ++trIt) {
		Tree * t = *trIt;
		if (t && NxsString::case_insensitive_equals(s.c_str(), t->GetName().c_str()))
			return t;
	}
	return NULL;
}

/*******************************************************************************
 * Adds the taxa in `newTax` to the current storage of tax labels. See the
 *	the comments on `UpdateMode` to see 
 */
const std::vector<std::string> & BullKernel::updateTaxa(const std::vector<std::string> & newTax, const UpdateMode mode)
	{
	openTaxa();
	const bool hadTaxa(!taxaLabels.empty());
	bool newTaxaAdded = false;
	std::map<std::string, unsigned> oldTaxLabelToInd;
	std::map<std::string, unsigned>::iterator msi;
	std::vector<std::string> scratchStringVec;
	if (mode == kClearUpdateMode)
		{
		taxaLabels.clear();
		taxLabelToInd.clear();
		if (hadTaxa)
			notifyListeners(BullListener::kTaxaDeleted);
		}
	else if (mode == kOnlyNewUpdateMode)
		{
		swap(oldTaxLabelToInd, taxLabelToInd);
		taxaLabels.clear();
		taxLabelToInd.clear();
		if (hadTaxa)
			notifyListeners(BullListener::kTaxaDeleted);
		for (std::vector<std::string>::const_iterator taxIt = newTax.begin(); taxIt != newTax.end(); ++taxIt)
			{
			std::string cap(NxsString::get_upper(*taxIt));
			if (oldTaxLabelToInd.find(cap) == oldTaxLabelToInd.end())
				{
				taxLabelToInd[cap] = taxaLabels.size();
				taxaLabels.push_back(*taxIt);
				newTaxaAdded = true;
				}
			}
		}
	else if (mode == kIntersectUpdateMode)
		{
		swap(oldTaxLabelToInd, taxLabelToInd);
		taxaLabels.clear();
		taxLabelToInd.clear();
		if (hadTaxa)
			notifyListeners(BullListener::kTaxaDeleted);
		for (std::vector<std::string>::const_iterator taxIt = newTax.begin(); taxIt != newTax.end(); ++taxIt)
			{
			const std::string cap(NxsString::get_upper(*taxIt));
			if (oldTaxLabelToInd.find(cap) != oldTaxLabelToInd.end())
				{
				taxLabelToInd[cap] = taxaLabels.size();
				taxaLabels.push_back(*taxIt);
				newTaxaAdded = true;
				}
			}
		}
	else if (mode == kReplaceUpdateMode)
		{
		taxaLabels.clear();
		taxLabelToInd.clear();
		if (hadTaxa)
			notifyListeners(BullListener::kTaxaDeleted);
		for (std::vector<std::string>::const_iterator taxIt = newTax.begin(); taxIt != newTax.end(); ++taxIt)
			{
			const std::string cap(NxsString::get_upper(*taxIt));
			taxLabelToInd[cap] = taxaLabels.size();
			taxaLabels.push_back(*taxIt);
			newTaxaAdded = true;
			}
		}
	else if (mode == kOnlyOldUpdateMode)
		{
		std::set<unsigned> indsToDel;
		for (std::vector<std::string>::const_iterator taxIt = newTax.begin(); taxIt != newTax.end(); ++taxIt)
			{
			const std::string cap(NxsString::get_upper(*taxIt));
			msi = taxLabelToInd.find(cap);
			if (msi != taxLabelToInd.end())
				indsToDel.insert(msi->second);
			}
		if (!indsToDel.empty())
			{
			taxLabelToInd.clear();
			scratchStringVec.swap(taxaLabels);
			taxaLabels.clear();
			for (unsigned i = 0; i < scratchStringVec.size(); ++i)
				{
				const std::string cap(NxsString::get_upper(scratchStringVec[i]));
				taxLabelToInd[cap] = taxaLabels.size();
				taxaLabels.push_back(scratchStringVec[i]);
				}
			notifyListeners(BullListener::kTaxaDeleted);
			}
		}
	else if (mode == kXorUpdateMode || mode == kAppendUpdateMode)
		{
		std::set<unsigned> indsToDel;
		std::vector<std::string> toAdd;
		for (std::vector<std::string>::const_iterator taxIt = newTax.begin(); taxIt != newTax.end(); ++taxIt)
			{
			const std::string cap(NxsString::get_upper(*taxIt));
			msi = taxLabelToInd.find(cap);
			if (msi == taxLabelToInd.end())
				toAdd.push_back(*taxIt);
			else if (mode == kXorUpdateMode)
				indsToDel.insert(msi->second);
			}
		if (!indsToDel.empty())
			{
			taxLabelToInd.clear();
			scratchStringVec.swap(taxaLabels);
			taxaLabels.clear();
			for (unsigned i = 0; i < scratchStringVec.size(); ++i)
				{
				const std::string cap(NxsString::get_upper(scratchStringVec[i]));
				taxLabelToInd[cap] = taxaLabels.size();
				taxaLabels.push_back(scratchStringVec[i]);
				}
			notifyListeners(BullListener::kTaxaDeleted);
			}
		for (std::vector<std::string>::const_iterator taxIt = toAdd.begin(); taxIt != toAdd.end(); ++taxIt)
			{
			const std::string cap(NxsString::get_upper(*taxIt));
			taxLabelToInd[cap] = taxaLabels.size();
			taxaLabels.push_back(*taxIt);
			newTaxaAdded = true;
			}
		}
	

	if (newTaxaAdded)
		notifyListeners(BullListener::kTaxaDeleted);
	closeTaxa();
	return taxaLabels;
	}


BullKernel::~BullKernel() {
	notifyListeners(BullListener::kDying);
	this->deleteModels(false);
	this->deleteTrees(false);
	notifyListeners(BullListener::kDead);
}

void BullKernel::Reset() {
	logEachStep = true; //TEMPORARY
	criterion = Criterion(MAX_LIKE);
	this->deleteModels(false);
	this->deleteTrees(false);
	workingOnBranches = false;
	analysisMode = AnalysisMode(INFERENCE);
	numSimulationConditionsDone = 0;
}


void BullKernel::deleteTrees(bool leaveOpen) {
	if (treeList.empty())
		return;
	openTrees();
	try {
		std::list<Tree *>::iterator trIt = this->treeList.begin();
		for (; trIt != this->treeList.end(); ++trIt)
			{
			Tree * trp = *trIt;
			delete trp;
			}
		this->treeList.clear();
	}
	catch (...) {
		if (!leaveOpen)
			closeTrees();
	}
}

void BullKernel::deleteModels(bool leaveOpen) {
	openModels();
	try {
		modelSet.freeMemory();
		}
	catch (...) {
		if (!leaveOpen)
			closeModels();
	}
}

////////////////////////////////////////////////////////////////////////////////
/// Assumes control of the trees in `treesToBeAdded` according to `updateMode`
/// \returns a vector of trees that the caller is responsible for deleting
////////////////////////////////////////////////////////////////////////////////
std::vector<Tree *> BullKernel::updateTrees(
  std::vector<Tree *> treesToBeAdded,
  const UpdateMode updateMode) { 
	if (updateMode == kIgnoreUpdateMode)
		return treesToBeAdded;

	BullAttributeManager bam(BullAttribute(kTrees), *this);

	if (updateMode == kClearUpdateMode || updateMode == kReplaceUpdateMode) {
		this->deleteTrees(true);
	}
	else if (updateMode != kAppendUpdateMode) {
		string msg = "Tree storage mode ";
		AppendInt(msg, (int) updateMode);
		msg.append(" is not currently supported.");
		throw XBull(msg);
	}
	
	if (updateMode == kClearUpdateMode)
		return treesToBeAdded;

	// Transfer ownership of the trees to treeList
	std::vector<Tree *> unused;
	for (std::vector<Tree *>::iterator i = treesToBeAdded.begin(); i != treesToBeAdded.end(); ++i)
		treeList.push_back(*i);

	notifyListeners(BullListener::kTaxaAdded);
	return unused;
}

void BullKernel::notifyListeners(BullListener::Event e, void *blob) {
	void * t = (void *) this;
	for (std::set<BullListener *>::iterator i = listeners.begin(); i != listeners.end(); ++i)
		(*i)->stateChanged(e, t, blob);
}

std::string bull::EncodingTypeToString(EncodingType et)
{
	if (et == EncodingType(DNANoGap))
		return std::string("dna");
	if (et == EncodingType(AminoAcid)) 
		return std::string("protein");
	assert(0);
	throw XBull("Unknown EncodingType");
}

bool BullKernel::hasModel() const 
{
	return !this->modelSet.modelPtrs.empty();
}
void BullKernel::simulate(const SimulateOpts & simOptions)
{	
 	analysisMode = AnalysisMode(SIMULATION);
	criterion = Criterion(MAX_LIKE);

	std::vector<EncodingType> outTypes = simOptions.outTypes;
	if (outTypes.empty())
		return;
	//TEMPORARY lots of the interaction with simSettings assumes that you are using SSRFCodonSubModel
	
	//alias the array of model pointers as modArr
	if (this->modelSet.modelPtrs.empty())
		throw XBull("You can't simulate a codon model without first declaring the models using codLikeStartVal");
	SSRFCodonSubMod **modArr = &(this->modelSet.modelPtrs[0]); 
	assert(modArr);
	assert(*modArr);
	const unsigned nCodons = (*modArr)->nCodonsInProtein;
	
	// set concatenations based on simOptions.concatenations or simOptions.nSimChars;
	unsigned concatenations = simOptions.concatenations;
	if (concatenations < 1) {
		assert (simOptions.nSimChars > 0);
		if ((simOptions.nSimChars)%(3*nCodons))
			throw XBull("The number of characters must be a multiple of the number of amino acids");
		concatenations = (simOptions.nSimChars)/(3*nCodons);
	}

	SimSettings simSettings;
	
	ModifyBlenMultSoBranchesAreScaled(modArr, nCodons);
	std::vector<Model *> subModelVec(1, 0L);
	for (unsigned k = 0; k < nCodons; k++) {
		subModelVec[0] = modArr[k];
		simSettings.AddPartSetting(PartSettings(k, k+1), subModelVec);
	}
	simSettings.setNumConcats(concatenations);
	simSettings.setNumReps(simOptions.nReps);
	const unsigned cnr = (concatenations * simOptions.nReps);
	const unsigned nsao = (cnr < 100 ? cnr : 100)/concatenations;
	simSettings.setNumSimsAtOnce(nsao);
	


	vector<std::string> pbfilename;
	simSettings.setNumOutputs(outTypes.size());
	for (unsigned i = 0; i < outTypes.size(); i++) {
		simSettings.SetOutputDatatype(i, outTypes[i]);
	}
	vector<string> outputFilenames = simOptions.outputFilenames;
	const bool automatic = (simOptions.automatic || outputFilenames.size() < outTypes.size());
	
	SimNodeLikeAttr::verboseSimulation =  bool(simOptions.verbose);
	

	vector<Tree *> modelTreesPtrVec;
	for (std::list<Tree *>::iterator trIt = treeList.begin(); trIt != treeList.end(); ++trIt)
		modelTreesPtrVec.push_back(*trIt);
	
	
	if (modelTreesPtrVec.empty())
		throw XBull("Model tree not loaded.");

	//Set up a file for the versions of the model tree with branches collapsed 
	// if there were no changes across the branch.
	ofstream collf;
	ofstream *collapseTreeFilePtr = NULL;
	ofstream aacollf;
	ofstream *aaCollapseTreeFilePtr = NULL;
	std::string aaCollapsedTreeFilename("aa_");
	if (!simOptions.collapsedTreeFilename.empty()) {
		aaCollapsedTreeFilename.append(simOptions.collapsedTreeFilename);
		collf.open(simOptions.collapsedTreeFilename.c_str());
		assert(collf.good());
		modelTreesPtrVec[0]->PrintTaxaBlock(collf);
		collf<<"endblock; \n\nbegin trees; \n";
		collapseTreeFilePtr = &collf;
		
		aacollf.open(aaCollapsedTreeFilename.c_str());
		assert(aacollf.good());
		modelTreesPtrVec[0]->PrintTaxaBlock(aacollf);
		aacollf << "endblock; \n\nbegin trees; \n";
		aaCollapseTreeFilePtr = &aacollf;
	}
				
	if (simOptions.seed != 0)
		bull::gRngSeed = simOptions.seed;
	for (std::size_t i = 0; i < modelTreesPtrVec.size(); i++) {
		//std::cout << "Tree " << i << std::endl;
		numSimulationConditionsDone++;
		Tree * modelTree = modelTreesPtrVec[i];
		assert(modelTree);
		if (!(modelTree->hasBranchLengths)) {
			//weak test (doesn't look at all partitions)
			throw XBull( "The tree must have branch lengths in order to simulate");
		}

		PrepareTreeForSimulation(*modelTree, simSettings);
		unsigned  ndone = 0;
		while (simOptions.nReps > ndone) {
			const unsigned ntodothisRound = std::min((unsigned)(simOptions.nReps-ndone), simSettings.getNumSimsAtOnce());
			long seedatbeg = bull::gRngSeed;
			modelTree->SimulateData(ntodothisRound); //done in batches to speed things up
			
			for (unsigned nofs = 0; nofs < outTypes.size(); nofs++) {
				std::string thisRepsOutFile;
				if (automatic) {
					thisRepsOutFile = EncodingTypeToString(outTypes[nofs]);
					if (i != 0)
						simSettings.setAppendToOut(true);
					if (!simOptions.tagname.empty())
						thisRepsOutFile += simOptions.tagname;
				}
				else {
					thisRepsOutFile = outputFilenames[nofs];
				}
				WriteSimulations(ntodothisRound, modelTree, thisRepsOutFile, outTypes[nofs], ndone, seedatbeg, numSimulationConditionsDone, simSettings);
			}
			if (collapseTreeFilePtr) {
				std::string oriname;
				oriname = modelTree->GetName();
				std::string temp = "CollapsedModel";
				modelTree->SetName(temp);
				modelTree->PrintBranchesWithMuts(*collapseTreeFilePtr, true, false);
				modelTree->PrintBranchesWithMuts(*aaCollapseTreeFilePtr, true, true);
				modelTree->SetName(oriname);
			}
			ndone += ntodothisRound;
		}

		DestroySimulationAttributes(*modelTree);
	}
	std::cout << "Simulated data written to: ";
	for (unsigned i = 0; i < outTypes.size(); i++) {
		std::cout << outputFilenames[i] << '\n';
	}
	if (collapseTreeFilePtr) {
		*collapseTreeFilePtr << "end; " << endl;
		*aaCollapseTreeFilePtr << "end; " << endl;
		std::cout << "Collapsed tree written to " << simOptions.collapsedTreeFilename << '\n';
		std::cout << "Amino Acid collapsed tree written to " << aaCollapsedTreeFilename << '\n';
	}
	collf.close();
	aacollf.close();
	
}

void BullKernel::WriteSimulations(
  unsigned num,
  Tree *t,
  std::string outpfn,
  unsigned outputEncodingType,
  unsigned totRepNum, 
  long seedAtBeg,
  unsigned simCondNum,
  const SimSettings & sSettings)
{
	ofstream outpf;
	if (!totRepNum && !(sSettings.getAppendToOut())) {
		outpf.open(outpfn.c_str());
		outpf<<"#NEXUS\n";
	}
	else {
		outpf.open(outpfn.c_str(), ios::app);
	}

	
	SetOfTreeSimAttr *tssa = t->GetSimAttributesPtr();
	unsigned totntax = t->GetNTax();
	unsigned nparts = tssa->GetNParts();
	unsigned totNCharPerDataSet = 0;
	vector<unsigned> ncharsInPart(nparts);

	//assumes that all characters will be stored in the simattributes as character state ordination codings in a short
	bool allsamemodelencoding = true;
	unsigned theModelEncodingType = tssa->GetModel(0, 0)->GetEncodingType();
	for (unsigned i = 1; i < nparts; i++) {
		if (tssa->GetModel(i, 0)->GetEncodingType()!= theModelEncodingType) {
			allsamemodelencoding = false;
			throw MTHException("Right now only one type of model encoding at a time can be used to simulate data");
		}
	}
	double modelToOutputCharRatio = NumDNACharactersPerCharacter(theModelEncodingType)*NumColumnsPerCharacter(outputEncodingType);
	modelToOutputCharRatio/= NumDNACharactersPerCharacter(outputEncodingType);
	for (unsigned i = 0; i < nparts; i++) {
		ncharsInPart[i] = tssa->GetNCharInEachDataSet(i);
		totNCharPerDataSet += (int)(((double) ncharsInPart[i])*modelToOutputCharRatio);
	}
	
	unsigned maxlen = t->GetMaxTaxonLabelLength(true)+4;
	char **OrdToNexTrans;
	OrdToNexTrans = GetOrdinationToNexusTranslator(theModelEncodingType, outputEncodingType);
	if (!totRepNum) {
		outpf<<"[!###SIMCONDITIONS = "<<simCondNum<<"##ModelTree = ";
		t->Print(outpf, true, false, false);
		outpf<<"]\n";
		outpf<<"BEGIN PAUP; \n\texecute preDataPaup.nex; \nEND; \n";
	}
	for (unsigned ds = 0; ds < num; ds++) {
		if(!ds)
			outpf<<"[!Seed At Beginning Of This batch of simulations = " <<seedAtBeg<<"]\n";
		outpf<<"[!######REP"<<totRepNum+ds<<"]\nBegin Data; \nDimensions NTAX = "<<totntax;
		outpf<<" NCHAR = "<<totNCharPerDataSet<<"; \n";
		outpf<<" FORMAT datatype = "<<NexusNameOfEncoding(outputEncodingType)<<" ; \n MATRIX\n";
		for (unsigned tax = 0; tax < t->GetNTax(); tax++) {
			Node *tempno = t->GetTerminalNode(tax);
			std::string tn = tempno->GetName();
			if (ds)
				outpf<<tn;
			else	
				outpf<<BlanksToUnderscores(tn);
			for (unsigned spa = tn.length(); spa < maxlen; spa++)
				outpf << " ";
			for (unsigned parti = 0; parti < nparts; parti++) {
				unsigned citp = ncharsInPart[parti];
				tempno->SetSimAtt(parti, 0); //TEMPORARY only works for one model per partition simulations
				const SimNodeLikeAttr * siminf = tempno->simInfo;
				assert(siminf);
				assert(siminf->charsAsShorts.size() >= (ds+1)*citp);
				const unsigned offset = ds*citp;
				for (unsigned ch = 0; ch < citp; ch++) {
					const unsigned globCodonNumber = siminf->getStateGlobalEncoding(offset + ch);
					outpf << OrdToNexTrans[ globCodonNumber ];
				}
			}
			outpf<<"\n";
		}
		outpf<<"; \nENDBLOCK; \n";
		outpf<<"BEGIN PAUP; \n\texecute eachDataPaup.nex; \nEND; \n";
		}

	// append a standard paup block with hooks

	outpf << "BEGIN PAUP; \n\texecute postDataPaup.nex; \nEND; \n";

	outpf.close();
}




/**
 * This frees the memory taken by a tree's sim attributes
 * after the simulation is done.  It is called by HandleSimulate
*/
void BullKernel::DestroySimulationAttributes(Tree &t)
{	
	Node ** const recursiveNodeList = t.GetRecursiveNodeList();
	assert(recursiveNodeList);
	Node **temp;
	SetOfTreeSimAttr *treesSetofSim = t.GetSimAttributesPtr();
	PartModIndex potOwner;
	if (treesSetofSim) {
		for (unsigned p = 0; p < treesSetofSim->GetNParts(); p++) {
			for (unsigned m = 0; m < treesSetofSim->GetNModels(p); m++) {
				potOwner.Set(p, m);
				if (!(treesSetofSim->GetOwnsBrLengths(potOwner))) {
					temp = recursiveNodeList;
					while (*temp) {
						(*temp)->SetSimAtt(potOwner.partNum, potOwner.modNum); //@mth this doesn't delete!
						(*temp)->simInfo->SetBLenPtr(NULL);
						temp++;
					}
				}
			}
		}
	}		
	temp = recursiveNodeList;
	while (*temp)
		{if((*temp)->IsTerm())
			(*temp)->ReplaceSetOfLikeInfoButKeepBranchLen(NULL); //Destroys old SimAttr for the nodes (including brLens)
		else
			(*temp)->ReplaceSetOfLikeInfoButKeepBranchLen(NULL); //Destroys old SimAttr for the nodes (including brLens)
		temp++;
		}
	t.SetSimAttributes(NULL); //Destroys old SimAttr for the tree (including brLen Modifiers)
}

/*
 * @method PrepareTreeForSimulation [void :protected]
 * @param t [Tree *] the tree to be scored
 * 
 * Replaces/instantiates like attributes for the tree and it nodes, 
 * initializes parameters, and if using the SSRF model EliminatesZeros 
 * called by BullShell::HandleSimulate()
 * 
*/
void BullKernel::PrepareTreeForSimulation(Tree &t, SimSettings & sSettings)
{	
	Node **temp = t.GetRecursiveNodeList();
	while (*temp)
		{
		(*temp)->ReplaceSetOfLikeInfoButKeepBranchLen(CreateSetOfSimAttrForANode(sSettings));
		temp++;
		}
	t.SetSimAttributes(CreateSetOfSimAttrForATree(sSettings));

}
/**
 * @method CreateSetOfSimAttrForANode [SetOfLikeAttr *:protected]
  *
 * DANGEROUS way to give an node its SimAttr and a setofLikeAttr for
 * the simulations.	 bool 
 * NOTE that SimAttr are derived from LikeAttr, and this returns a 
 * MEMORY ALLOCATION for like attr and setOfLikeAttr.  Deletion occurs in
 * Node::ReplaceSetOfLikeInfoButKeepBranchLen( SetOfLikeAttr *) which is called by BullShell::PrepareTreeForSimulation()
 * with CreateSetOfSimAttrForANode as an argument
*/
SetOfLikeAttr *BullKernel::CreateSetOfSimAttrForANode(SimSettings & sSettings)
{	//For simulations all nodes get a Terminal Like attr (because the characters will be stored as shorts
	//if memory is an issue the internal nodes could share memory
	SetOfLikeAttr *temp;
	temp = new SetOfLikeAttr(sSettings.GetNDataPartitions());
	const unsigned ncnsao = sSettings.getNumConcats() * sSettings.getNumSimsAtOnce();
	for (unsigned i = 0; i < sSettings.GetNDataPartitions(); i++) {
		temp->SetNumModsInPart(i, 1);
		Model * tempMod = sSettings.GetModel(i, 0);
		unsigned len =  ncnsao*sSettings.GetPartitionLength(i);
		SimNodeLikeAttr * snla = new SimNodeLikeAttr(len, tempMod);
		temp->AddLikeAtt(i, 0, snla);
	}
	return temp;
}

/**
 * klutzy way to give an tree its SimAttr and a setofLikeAttr for
 * the simulations. 
 * MEMORY ALLOCATION for like attr and setOfLikeAttr.  Deletion occurs in
 * Tree::SetSimAttributes(SetOfTreeSimAttr *) which is called by BullShell::PrepareTreeForSimulation()
 * with CreateSetOfSimAttrForATree as an argument
*/
SetOfTreeSimAttr *BullKernel::CreateSetOfSimAttrForATree(SimSettings & sSettings)
{	
	SetOfTreeSimAttr * temp = new SetOfTreeSimAttr(sSettings.GetNDataPartitions());
	for (unsigned i = 0; i < sSettings.GetNDataPartitions(); i++) {
		temp->SetNumModsInPart(i, 1);
		const unsigned singleLen = sSettings.getNumConcats() * sSettings.GetPartitionLength(i);
		TreeSimAttr *tempPtr;
		const unsigned len = singleLen * sSettings.getNumSimsAtOnce();
		tempPtr = new TreeSimAttr(len, *sSettings.GetModel(i, 0));
		temp->AddAtt(i, 0, tempPtr);
#		ifdef ONESETOFBRANCHES
			tempPtr->SetOwnsBrLen(i == 0);
#		else
			MTHException("Unwritten code in CreateSetOfSimAttrForATree()");
#		endif
		temp->SetNCharInEachDataSet(i, singleLen);
	}
	return temp;
}			



