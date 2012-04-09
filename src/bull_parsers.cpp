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

#include "bull_parsers.hpp"
#include <sstream>
#include <algorithm>

#include "bull.hpp"
#include "xbull.hpp"
#include "string_extensions.hpp"
#include "ncl/nxstreesblock.h"

using namespace std;
using namespace bull;

int GetHighestFileNum(std::string tag,int num=0);
void ReadNexusVector(std::vector<std::string> *v, NexusToken &token);

void SkipEqualsIfNext(NexusToken & token);

int GetHighestFileNum(std::string tag,int num/*=0*/) {
	ifstream testf;
	std::string tempstr;
	tempstr=tag;
	tempstr+=num;
	testf.open(tempstr.c_str());
	while (testf.good()) {
		testf.close();
		num++;
		tempstr=tag;
		tempstr+=num;
		testf.open(tempstr.c_str());
	}
	testf.close();
	return num - 1;
}



void BullShell::HandleGetTrees(NexusToken& token ) {
	GetTreesOpts gto = ParseGetTrees(token);
	try {
		ExecGetTrees(gto);
	}
	catch (XBull & x) {
		throw XBull(x.msg, token);
	}
}

void BullShell::HandleCodLikeStartVal(NexusToken& token ) {
	CodLikeStartOpts clso;
	ParseCodLikeStartValCommand(clso, token);
	try {
		ExecCodLikeStartValCommand(clso);
	}
	catch (XBull & x) {
		throw XBull(x.msg, token);
	}
}

void BullShell::ExecCodLikeStartValCommand(const CodLikeStartOpts &clso)
{
	SSRFCodonSubModSet s;
	s.initialize(clso.gtrParams, clso.aaFreqs, clso.multipliers, clso.treeScale, clso.genetic_code);
	kernel.setModelConstRef(s);
	s.surrenderThenClear();
	assert(kernel.hasModel());
}

void BullShell::HandleSimulate(NexusToken& token ) {
	SimulateOpts sto = ParseSimulateCommand(token);
	try {
		ExecSimulateCommand(sto);
	}
	catch (XBull & x) {
		throw XBull(x.msg, token);
	}
}

/// Throws an XBull excpetion if token is not ";"
void BullShell::RequireSemicolon(NexusToken & token, const char *cmd) const {
	if ( !token.Equals(";"))
		ThrowNoSemicolon(token, cmd);
}

void BullShell::ThrowNoSemicolon(NexusToken & token, const char *cmd) const {
	errormsg = "Expecting ';'";
	if (cmd != NULL) {
		errormsg += " to terminate the ";
		errormsg.append(cmd);
		errormsg += "command,";
	}
	errormsg += " but found ";
	errormsg += token.GetToken();
	errormsg += " instead";
	throw XBull(errormsg, token);
}
/**
 * @method HandleEndblock [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XBull
 *
 * Called when the END or ENDBLOCK command needs to be parsed
 * from within the BullShell block.	Basically just checks to make
 * sure the next token in  the data file is a semicolon.
 */
void BullShell::HandleEndblock( NexusToken& token ) {
	token.GetNextToken();
	RequireSemicolon(token, "END");
}

/**
 * @method HandleExecute [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XBull
 *
 * Handles everything after the EXecute keyword and the terminating
 * semicolon.  Flushes all blocks before executing file specified,
 * and no warning is given of this
 */
void BullShell::HandleExecute( NexusToken& token ) {
	token.GetNextToken();
	std::string fn = token.GetToken();
	token.GetNextToken();
	RequireSemicolon(token, "EXECUTE");
	cout << "\nExecuting " << fn << "..." << endl;
	readFile(fn, false);
	cout << "Done with file "<<fn<<endl;
}

GetTreesOpts BullShell::ParseGetTrees(NexusToken& token ) const {
	errormsg.clear();
	if (kernel.getNumTaxa() == 0) {
		errormsg="You can't get trees without an active set of taxa";
		throw XBull( errormsg, token);
	}
	GetTreesOpts gto;
	bool readingMostRecent = false;
	gto.fromTree = -1;
	gto.toTree = -1;
	gto.mode = 7;
	std::string filePref;

	token.GetNextToken();
	while (token.GetToken()!=";") {
		if (token.Abbreviation("STOREBrlens")) {
			token.GetNextToken();
			gto.storeBrLensFromFile=true;
			if (token.GetToken() == "=") {
				token.GetNextToken();
				if (token.Abbreviation("No"))
					gto.storeBrLensFromFile=false;
				else if (!token.Abbreviation("Yes")) {
					errormsg="Expecting YES or NO after CurrentBranchLengths = option to LSCORE command";
					throw XBull(errormsg, token);
				}
			}
		}
		else if (token.Abbreviation("FRom")) {
			DemandEquals(token, "after From subcommand of GetTrees command");
			gto.fromTree = DemandPositiveInt(token, "after From subcommand of GetTrees command");
			if (gto.fromTree < 0) {
				errormsg="Expecting positive integer after FROM option of GETTREES";
				throw XBull(errormsg, token);
			}	
		}
		else if (token.Abbreviation("To")) {
			DemandEquals(token, "after To subcommand of GetTrees command");
			gto.toTree = DemandPositiveInt(token, "after To subcommand of GetTrees command");
			if (gto.toTree < 0) {
				errormsg="Expecting positive integer after TO option of GETTREES";
				throw XBull(errormsg, token);
			}	
		}
		else if (token.Abbreviation("REplace"))
			gto.mode=3;
		else if (token.Abbreviation("Mode")) {
			DemandEquals(token, "after Mode subcommand of GetTrees command");
			gto.mode = DemandPositiveInt(token, "after Mode subcommand of GetTrees command");
			if (gto.mode!=3 && gto.mode!=7) {
				errormsg="Right now bull only get trees with mode 3 or mode 7";
				throw XBull(errormsg, token);
			}	
		}
		else if (token.Abbreviation("FIle")) {
			DemandEquals(token, "after File subcommand of GetTrees command");
			token.GetNextToken();
			gto.filename = token.GetTokenReference();
		}
		else if (token.Abbreviation("PREfix")) {
			token.GetNextToken();
			if (token.GetToken() == "=")	
				token.GetNextToken();
			filePref = token.GetToken();
		}
		else if (token.Abbreviation("MOStrecent"))
			readingMostRecent=true;
		else {
			errormsg = "Unrecognized subcommand (";
			errormsg << token.GetTokenReference() << ") in GetTrees command.";
			throw XBull(errormsg, token);
		}
		token.GetNextToken();
	}
	
	if (readingMostRecent) {
		int highNum = GetHighestFileNum(filePref);
		if (highNum < 0) {
			errormsg << "Couldn't open file " << filePref << "0";
			throw XBull(errormsg, token);
		}  
		gto.filename = filePref;
		gto.filename += highNum;
	}
	if (gto.fromTree > -1 && gto.toTree > gto.fromTree) {
		errormsg << "The FromTree setting (" << gto.fromTree << ") cannot be less than the ToTree setting (" << gto.toTree << ")";
		throw XBull(errormsg, token);
	}
		
	gto.fromTree = (gto.fromTree < 1 ? -1 : gto.fromTree - 1);
	gto.toTree = (gto.toTree < 1 ? -1 : gto.toTree - 1);
	return gto;
}

void BullShell::ExecGetTrees(const GetTreesOpts & gto) {
	errormsg.clear();
	BlockReaderList blocks = getBlocksFromFile(gto.filename, true);
	if (blocks.empty()) {
		errormsg << "No Trees Block found in " << gto.filename;
		throw XBull(errormsg);
	}
	
	unsigned ntreesInFile = 0;
	const unsigned ntreesBlocks = blocks.size();
	vector<NxsTreesBlockAPI *> trBlocks;
	for (BlockReaderList::iterator b = blocks.begin(); b != blocks.end(); ++b) {
		NxsTreesBlockAPI *tb = (NxsTreesBlockAPI *)(*b);
		trBlocks.push_back(tb);
		assert(tb);
		ntreesInFile += tb->GetNumTrees();
	}

	if (ntreesInFile < 1) {
		errormsg << "No Trees found in " << gto.filename;
		throw XBull(errormsg);
	}
	
	const unsigned toTree = (unsigned) (gto.toTree < 0 ? 0 : gto.toTree); 
	const unsigned fromTree = (unsigned) (gto.fromTree < 0 ? ntreesInFile - 1 : gto.fromTree);
	if (fromTree > toTree) {
		errormsg << "Cant get Trees from " << fromTree << " to " << toTree;
		throw XBull( errormsg);
	}

	// Store the specified trees
	// these trees are deleted in the for loop if there is an error, otherwise
	//	they are put under the control of the kernel
	vector<Tree *>	treesToBeAdded; 
	vector<NxsTreesBlockAPI *>::const_iterator tbIt = trBlocks.begin();
	NxsTreesBlockAPI * currTB = *tbIt;
	unsigned nTreesInThisBlock = currTB->GetNumTrees();
	int nleftInThisBlock = (int) nTreesInThisBlock;
	for (unsigned i = 0; i < toTree; i++, --nleftInThisBlock) {
		while (nleftInThisBlock < 1) {
			++tbIt;
			assert(tbIt != trBlocks.end());
			currTB = *tbIt;
			nTreesInThisBlock = currTB->GetNumTrees();
			nleftInThisBlock = (int) nTreesInThisBlock;
		}
		if (i >= fromTree) {
			const std::string newick = currTB->GetTranslatedTreeDescription((unsigned)(nTreesInThisBlock - nleftInThisBlock));
			Tree * temptree = new Tree(newick, gto.storeBrLensFromFile);
			if (!temptree->IsGood()) {
				for (unsigned j = 0; j < i - fromTree; j++)
					delete treesToBeAdded[j];
				errormsg << "Problem Reading Tree Description of " << currTB->GetTreeName(i);
				delete temptree;
				throw XBull(errormsg);
			}
			std::string s = currTB->GetTreeName(i);
			ToUpper(s);
			temptree->SetName(s);
			treesToBeAdded.push_back(temptree);
		}
	}

	ioObject.message.clear();
	ioObject.message << (int)(toTree - fromTree) << " trees read from " << ntreesBlocks << " TREES block(s) in " <<  gto.filename;
	ioObject.printMessage(BullIO::STATUS_MSG_LEVEL);

	assert(gto.mode <= 7 && gto.mode >= 0);
	kernel.updateTrees(treesToBeAdded, UpdateMode(gto.mode));
	for (BlockReaderList::iterator b = blocks.begin(); b != blocks.end(); ++b) {
		if (*b != this)
			delete *b;
	}
}		


/**
 * @method HandleLog [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XBull
 *
 * Called when the LOG command needs to be parsed
 * from within the BullShell block.
 */
void BullShell::HandleLog( NexusToken& token )
{
	bool starting = false;
	bool stopping = false;
	bool appending = false;
	bool replacing = false;
	std::string logfname;
	errormsg.clear();
	// Retrieve all tokens for this command, stopping only in the event
	// of a semicolon or an unrecognized keyword
	//
	for (token.GetNextToken(); !token.Equals(";"); token.GetNextToken()) {
		if ( token.Abbreviation("STOp") )
			stopping = true;
		else if ( token.Abbreviation("STArt") )
			starting = true;
		else if ( token.Abbreviation("Replace") )
			replacing = true;
		else if ( token.Abbreviation("Append") )
			appending = true;
		else if ( token.Abbreviation("File") ) {
			DemandEquals(token, "after file subcommand of Log command");
			token.GetNextToken();
			logfname = token.GetTokenReference();
		}
		else {
			errormsg << "Unexpected keyword (" << token.GetToken() << ") encountered reading LOG command";
			throw XBull( errormsg, token);
		}
   }
   
   // Check for incompatible combinations of keywords
   if ( stopping && ( starting || appending || replacing || !logfname.empty() ) ) {
		errormsg = "Cannot specify STOP with any of the following START, APPEND, REPLACE, FILE";
		throw XBull( errormsg, token);
   }
   
   if ( appending && replacing ) {
		errormsg = "Cannot specify APPEND and REPLACE at the same time";
		throw XBull( errormsg, token);
   }

   if ( stopping ) {
		ioObject.stopLogging();
		kernel.logEachStep=false;
		return;
   }
   
	if (logfname.empty() ) {
		errormsg << "Must provide a file name when opening a log file\ne.g., log file=doofus.txt start replace;";
		throw XBull(errormsg, token);
	}	
	if (!appending && !replacing && FileExists(logfname.c_str())) {
		errormsg << "File " << logfname << " exists.\n APPEND or REPLACE must be specified in the LOG command.";
		throw XBull( errormsg, token);
	}
	ioObject.startLogging(logfname,  appending);	
}



void ReadNexusVector(std::vector<std::string> *v, NexusToken &token) { // TAH changed from above
	if (token.Equals("(")) {
		token.GetNextToken();
		while (!token.Equals(")")) {
			v->push_back(token.GetToken());
			token.GetNextToken();
		}
	}
	else 
		v->push_back(token.GetToken());
}


EncodingType BullShell::InterpretOutputType(NexusToken &token) const {
	if (token.Abbreviation("PROtein") || token.Abbreviation("AMino"))
		return EncodingType(AminoAcid);
	if (token.Abbreviation("Dna"))
		return EncodingType(DNANoGap);
	errormsg.clear();
	errormsg << "The output type " << token.GetToken() << " is invalid";
	throw XBull(errormsg, token);
}


SimulateOpts BullShell::ParseSimulateCommand(NexusToken& token ) const {	
	token.GetNextToken();
	int nOutputs = 0;
	int nSimChars = 0;
	SimulateOpts simOpts;
	const unsigned nTrees = kernel.getNumTrees();
	for (; token.GetToken() != ";"; token.GetNextToken()) { 
		if (token.Abbreviation("NReps")) {
			DemandEquals(token, "after NReps subcommand of Simulate command");
			simOpts.nReps = DemandPositiveInt(token, "after NReps subcommand of Simulate command");
		}
		if (token.Abbreviation("Seed")) {
			DemandEquals(token, "after Seed subcommand of Simulate command");
			simOpts.seed = DemandPositiveInt(token, "after Seed subcommand of Simulate command");
		}
		else if (token.Abbreviation("COLLapsefile")) {
			DemandEquals(token, "after CollapseFile subcommand of Simulate command");
			token.GetNextToken();
			simOpts.collapsedTreeFilename = token.GetTokenReference();
		}
		else if (token.Abbreviation("NChars")) {
			DemandEquals(token, "after NChars subcommand of Simulate command");
			simOpts.nSimChars = DemandPositiveInt(token, "after NChars subcommand of Simulate command");
		}
		else if (token.Abbreviation("Concat")) {
			DemandEquals(token, "after Concat subcommand of GetTrees command");
			simOpts.concatenations  = DemandPositiveInt(token, "after Concat subcommand of GetTrees command");
		}
		else if (token.Abbreviation("Verbose")) {
			simOpts.verbose = true;
		}
		else if (token.Abbreviation("NOUTput")) {
			DemandEquals(token, "after NOutput subcommand of Simulate command");
			nOutputs = DemandPositiveInt(token, "after NOutput subcommand of Simulate command");
			if (nOutputs < 1) {
				errormsg="nOutputs must be > 0";
				throw XBull( errormsg, token);
			}
		}
		else if (token.Abbreviation("Paupblockfile")) {
			DemandEquals(token, "after Paupblockfile subcommand of Simulate");
			token.GetNextToken();
			ReadNexusVector(&(simOpts.paupBlockFilename), token);
		}
		else if (token.Abbreviation("OVERWrite"))
			simOpts.overwrite = true;
		else if (token.Abbreviation("File")) {
			DemandEquals(token, "after File subcommand of Simulate");
			token.GetNextToken();
			ReadNexusVector(&(simOpts.outputFilenames), token);
		}
		else if (token.Abbreviation("TAg")) {
			DemandEquals(token, "after Tag subcommand of Simulate command");
			token.GetNextToken();
			simOpts.tagname = token.GetTokenReference();
		}
		else if (token.Abbreviation("AUTOmatic"))
			simOpts.automatic = true;
		else if (token.Abbreviation("OUTputtypes")) {
			DemandEquals(token, "after OutputTypes subcommand of Simulate");
			token.GetNextToken();
			EncodingType et;
			if (token.Equals("(")) {
				token.GetNextToken();
				while (!token.Equals(")")) {
					et = InterpretOutputType(token);
					simOpts.outTypes.push_back(et);
					token.GetNextToken();
				}
			}
			else {
				et = InterpretOutputType(token);
				simOpts.outTypes.push_back(et);
			}
		}
		else {
			Tree * temptree;
			long tn;
			if (NxsString::to_long(token.GetTokenAsCStr(), &tn)) {
				if ((unsigned) tn > nTrees || tn < 1) {
					errormsg = "Tree ";
					errormsg += token.GetToken();
					errormsg += " unknown.	simulation aborted.";
					throw XBull( errormsg, token);
				}
				temptree = kernel.getTree(tn - 1);
			}
			else {
				try {
					temptree = kernel.FindTreeFromName(token.GetToken());
				}
				catch (NoSuchTree) {
					errormsg << "Tree " << token.GetToken() << " unknown.";
					throw XBull( errormsg, token);
				}
			}
			simOpts.treeNames.push_back(token.GetToken()); //keeps track of what the user called the tree
			simOpts.treeAlias.push_back(temptree);
		}	
	}


	if (simOpts.outTypes.size() == 0){
		simOpts.outTypes.push_back(EncodingType(DNANoGap));
	}
	if ((nOutputs > 0) && simOpts.outTypes.size() != (unsigned) nOutputs) {
		errormsg = "The number of OutputTypes in a Simulate command must match the number of specified in the NOutput option";
		throw XBull( errormsg, token);
	}
	if (simOpts.automatic) {
		if (simOpts.tagname.empty()) {
			errormsg = "The Tag option of the Simulate command must be used if the Automatic naming mode is requested.";
			throw XBull( errormsg, token);
		}
		if (!simOpts.outputFilenames.empty()) {
			errormsg = "The FILE option of the Simulate command cannot be used if the Automatic naming mode is requested (use the Tag instead of file).";
			throw XBull( errormsg, token);
		}
		simOpts.outputFilenames.clear();
		std::vector<EncodingType>::const_iterator ot = simOpts.outTypes.begin();
		for (; ot != simOpts.outTypes.end(); ot++) {
			string ofname = EncodingTypeToString(*ot);
			ofname.append(simOpts.tagname);
			simOpts.outputFilenames.push_back(ofname);
		}
	}
	else if (simOpts.outTypes.size() != simOpts.outputFilenames.size()) {
		errormsg = "The number of OutputTypes in a Simulate command must match the number of Files specified";
		throw XBull( errormsg, token);
	}
	if ((!simOpts.paupBlockFilename.empty()) && simOpts.paupBlockFilename.size() != simOpts.outTypes.size()) {
		errormsg = "The number of PaupBlockFile in a Simulate command must match the number of OutputTypes specified";
		throw XBull( errormsg, token);
	}
	

	//TEMPORARY lots of the interaction with SSettings assumes that you are using SSRFCodonSubModel
	if (!kernel.hasModel()) {
		errormsg="You can't simulate a codon model without first declaring the models using codLikeStartVal";
		throw XBull( errormsg, token);
	}
	
	if ( simOpts.concatenations == 0) {
		if (nSimChars == 0) {
			errormsg="Either the number of characters or the number of concatenations must be specified";
			throw XBull( errormsg, token);
		}
		const unsigned naa = kernel.getModelConstRef().getNumAASites();
		const unsigned ndna = 3*naa;
		if (nSimChars % ndna) {
			errormsg="The number of characters must be a multiple of the number of amino acids";
			throw XBull( errormsg, token);
		}
		simOpts.concatenations = ndna/nSimChars;
	}
	return simOpts;
}


void SkipEqualsIfNext(NexusToken & token)
{
	token.GetNextToken();
	if (token.Equals("="))
		token.GetNextToken();
}

void BullShell::ParseCodLikeStartValCommand(CodLikeStartOpts &opts, NexusToken& token) const
{
	opts.aaFreqs.clear();
	DblVector blankRow(21, 1.0/20.0);
	blankRow[20] = 0.0; //zero out the freq of the stop codon.
			
	opts.genetic_code = GenCode(MITO);
	token.GetNextToken();
	opts.gtrParams.assign(N_MUT_PARAMS, .25);
	while (token.GetToken()!=";") {
		if (token.Abbreviation("CURREntbranchlengths") || token.Abbreviation("CURRBranchlengths")) {
			token.GetNextToken();
			opts.currbrlen = true;
			if (token.Equals("=")) {
				token.GetNextToken();
				if (token.Abbreviation("No"))
					opts.currbrlen = false;
				else if (!token.Abbreviation("Yes"))
					throw XBull("Expecting YES or NO after CurrentBranchLengths = option to CodLikeStartVal command");
				token.GetNextToken();
			}
		}
		else if (token.Abbreviation("GEneticCode")) {
			SkipEqualsIfNext(token);
			if (token.Abbreviation("Mitochondrial"))
				opts.genetic_code = GenCode(MITO); 
			else if (token.Abbreviation("Nuclear"))
				opts.genetic_code = GenCode(NUCLEAR);
			else 
				throw XBull("Expecting either Mito or Nuclear code");
			}
		else if (token.Abbreviation("BASEfreq")) {
			SkipEqualsIfNext(token);
			if (!token.Equals("("))
				throw XBull("Expecting ( after basefreq option to CodLikeStartVal command");
			for (unsigned i = 0; i < 3; i++) {
				token.GetNextToken();
				opts.gtrParams[i] = atof(token.GetToken().c_str());
				if (opts.gtrParams[i] <= 0.0)
					throw XBull("basefreqs must be > 0.0 CodLikeStartVal command");
			}
			if (opts.gtrParams[1] + opts.gtrParams[2] + opts.gtrParams[0] >= 1.0)
				throw XBull("Sum of A, C, and G must be <1.0 CodLikeStartVal command");
			token.GetNextToken();
			if (!token.Equals(")"))
				throw XBull("Expecting ) after basefreq option to CodLikeStartVal command");
		}
		else if (token.Abbreviation("RMATrix")) {
			SkipEqualsIfNext(token);
			if (!token.Equals("("))
				throw XBull("Expecting ( after rmatrix option to CodLikeStartVal command");
			for (int i = 0; i < 6; i++) {
				token.GetNextToken();
				opts.gtrParams[3 + i] = atof(token.GetToken().c_str());
				if (opts.gtrParams[3 + i] <= 0.0) 
					throw XBull("rmatrix  must be >0.0 CodLikeStartVal command");
			}
			token.GetNextToken();
			if (!token.Equals(")")) 
				throw XBull("Expecting ) after rmatrix	 option to CodLikeStartVal command");
			}
#		ifdef	ALLOWMULTIHITS
			else if (token.Abbreviation("DOublehit")) { //not a great name, but allows multiplier abbreviation to work
				SkipEqualsIfNext(token);
				opts.gtrParams[MULTI_HIT_PARAM_INDEX] = atof(token.GetToken().c_str());
				if (opts.gtrParams[MULTI_HIT_PARAM_INDEX] < 0.0)
					throw XBull("DOublehit (MultipleHitProb)  must be >= 0.0 CodLikeStartVal command");
			}
#		endif
		else if (token.Abbreviation("AAFreq")) {
			SkipEqualsIfNext(token);
			if (!token.Equals("("))
				throw XBull("Expecting ( after aafreq option to CodLikeStartVal command");
			token.GetNextToken();
			long tmp;
			if (!NxsString::to_long(token.GetTokenAsCStr(), &tmp))
				throw XBull("Expecting number of codons after aafreq=( option to CodLikeStartVal command");
			if (tmp < 1)
				throw XBull("number of codons must be >0 in	 CodLikeStartVal command");

			unsigned naa = (unsigned) tmp;
			opts.aaFreqs.assign(naa, blankRow);
			for (unsigned ii = 0; ii < naa; ii++) {
				token.GetNextToken();
				if (!token.Equals("("))
					throw XBull("Expecting ( for next site in aafreq option to CodLikeStartVal command", token);
				DblVector & row = opts.aaFreqs[ii];
				for (int j=0; j < 20; j++) {
					token.SetLabileFlagBit(NxsToken::hyphenNotPunctuation);
					token.GetNextToken();
					row[j] = atof(token.GetToken().c_str());
					if (row[j]<0.0)
						throw XBull("each amino acid freq	must be >0.0 CodLikeStartVal command", token);
				}
				token.GetNextToken();
				if (!token.Equals(")")) 
					throw XBull("Expecting ) for next site in aafreq option to CodLikeStartVal command", token);
			}
			token.GetNextToken();
			if (!token.Equals(")"))
				throw XBull("Expecting ) after aafreq	option to CodLikeStartVal command");
		}
		else if (token.Abbreviation("MUlt")) {
			SkipEqualsIfNext(token);
			if (!token.Equals("(")) 
				throw XBull("Expecting ( after mult option to CodLikeStartVal command");
			if (opts.aaFreqs.size() < 1)
				throw XBull("number of codons must be specified (in the AAFreq subcommand) before the mult option is used in CodLikeStartVal command");
			opts.multipliers.assign(opts.aaFreqs.size(), 1.0);
			for (unsigned ii = 0; ii < opts.aaFreqs.size(); ii++) {
				token.GetNextToken();
				opts.multipliers[ii] = atof(token.GetToken().c_str());
				if (opts.multipliers[ii] < 0.0)
					throw XBull("multiplier	 must be > 0.0 CodLikeStartVal command");
			}
			token.GetNextToken();
			if (!token.Equals(")")) 
				throw XBull("Expecting ) after mult	 option to CodLikeStartVal command");
		}
		else if (token.Abbreviation("TReescale")) {
			SkipEqualsIfNext(token);
			opts.treeScale = atof(token.GetToken().c_str());
			if (opts.treeScale < 0.0)
				throw XBull("The Tree scaling	must be >0.0 CodLikeStartVal command");
		}
		token.GetNextToken();
	}

}
