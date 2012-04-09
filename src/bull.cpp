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
#include "bull.hpp"

#include <sstream>
#include <fstream>
#include <algorithm>

#include "ncl/nxspublicblocks.h"

#include "bull_parsers.hpp"
#include "string_extensions.hpp"
#include "tree.hpp"
#include "xbull.hpp"

extern long bull::gRngSeed;

double DistanceBasesMoved(double *curr,double *prev);
double DistanceRateParamsMoved(double *curr,double *prev);
double EuclideanDistance(int arrSize,double *f,double *s);
void NormalizeGTRParams(double *curr,double *currNorm);

const double bull::DEFAULT_START_EDGE_LEN = 0.1;
const double bull::BULL_SMALL_DBL = 1.0e-10;
const double bull::BULL_BIG_DBL = 1.0e303;


using namespace std;
using namespace bull;


BullShell::BullShell(BullKernel & kern)
  :kernel(kern) 
{
	id = "BULL";
}

void BullShell::processAssumptionsBlock(NxsAssumptionsBlock *) 
{
}

void BullShell::processCharactersBlock(NxsCharactersBlockAPI *) 
{
}

void BullShell::processTaxaBlock(NxsTaxaBlockAPI *taxa)
{
	unsigned ntaxa = taxa->GetNumTaxonLabels();
	if (ntaxa > 1)
		cout << ntaxa <<" taxa in from taxa block"<<endl;
	std::vector<std::string> labels = taxa->GetAllLabels();
	std::vector<std::string> allL = kernel.updateTaxa(labels);
	
	ioObject.message << ntaxa <<" taxa in from taxa block.\n";
	ioObject.message << allL.size() <<" tree(s) currently in memory.\n";
	ioObject.printMessage(BullIO::INFO_MSG_LEVEL);
}

void BullShell::processTreesBlock(NxsTreesBlockAPI *trees) 
{
	unsigned ntrees=trees->GetNumTrees();
	std::vector<Tree *> vt(ntrees, 0L);
	try
		{
		for (unsigned i = 0; i < ntrees; i++) {
			NxsString newick = trees->GetTranslatedTreeDescription(i);
			Tree *temptree= new Tree(newick);
			if (!temptree->IsGood()) {
				errormsg.clear();
				errormsg << "Failed to parse tree number " << i+1<< " in the trees block.";
				throw MTHException(errormsg.c_str());
			}
			std::string s = trees->GetTreeName(i);
			temptree->SetName(s);
			vt[i] = temptree;
			}
		}
	catch (...)
		{
		for (std::vector<Tree *>::iterator v = vt.begin(); v != vt.end(); ++v)
			{
			Tree *t = *v;
			if (t)
				delete t;
			else
				break;
			}
		throw;
		}
	std::vector<Tree *> allt = kernel.updateTrees(vt);

	ioObject.message << ntrees <<" tree(s) read from file.\n";
	ioObject.message << allt.size() <<" tree(s) currently in memory.\n";
	ioObject.printMessage(BullIO::INFO_MSG_LEVEL);
}

void BullShell::processNexusBlocks(BlockReaderList &blocks)
{
	for (BlockReaderList::iterator uIt = blocks.begin(); uIt != blocks.end(); ++uIt) {
		NxsBlock * b = * uIt;
		if (b)
			{
			NxsString bId = b->GetID();
			if (bId.EqualsCaseInsensitive("CHARACTERS") || bId.EqualsCaseInsensitive("DATA"))
				{
				NxsCharactersBlockAPI * ncb	 = static_cast<NxsCharactersBlockAPI *>(b);
				processCharactersBlock(ncb);
				}
			else if (bId.EqualsCaseInsensitive("TAXA"))
				{
				NxsTaxaBlockAPI * taxb	= static_cast<NxsTaxaBlockAPI*>(b);
				processTaxaBlock(taxb);
				}
			else if (bId.EqualsCaseInsensitive("TREES"))
				{
				NxsTreesBlockAPI * trb	= static_cast<NxsTreesBlockAPI*>(b);
				processTreesBlock(trb);
				}
			else if (bId.EqualsCaseInsensitive("ASSUMPTIONS") 
					|| bId.EqualsCaseInsensitive("SETS")
					|| bId.EqualsCaseInsensitive("CODONS"))
				{
				NxsAssumptionsBlock * assump  = static_cast<NxsAssumptionsBlock*>(b);
				processAssumptionsBlock(assump);
				}
			}
	}
}

void BullShell::identifyAndProcessUnhandledBlocks(BlockReaderList &blocks)
{
	BlockReaderList unhandled;
	for (BlockReaderList::reverse_iterator uIt = blocks.rbegin(); uIt != blocks.rend(); ++uIt) {
		if (*uIt == this)
			break;
		unhandled.push_front(*uIt);
	}
	if (!unhandled.empty())
		this->processNexusBlocks(unhandled);
}

/*******************************************************************************
 * Checks the stack of NEXUS readers for public blocks that have been read since
 *	the last BULL block and then calls identifyAndProcessUnhandledBlocks()
 *	to handle them.
 */
void BullShell::processInterveningPublicNexusBlocks()
{
	if (readerStack.empty())
		return;
	NxsReader * reader = readerStack.top();
	if (!reader)
		return;
	BlockReaderList blocks = reader->GetUsedBlocksInOrder();
	identifyAndProcessUnhandledBlocks(blocks);
}

BlockReaderList BullShell::getBlocksFromFile(const std::string &filepath, const bool treesBlockOnly) 
{
	if (filepath.empty())
		throw NxsException ("Invalid (empty) filename specified", 0, 0, 0);
	ifstream inf(filepath.c_str(), ios::binary);
	if (!inf.good()) {
		NxsString err;
		err << "Could not read the file \"" << filepath <<"\"";
		throw NxsException(err, 0, 0, 0);
	}
	NxsToken token(inf);	

	ExceptionRaisingNxsReader nexusReader;
	NxsCloneBlockFactory factory;
	nexusReader.AddFactory(&factory);

	NxsAssumptionsBlock assumpB(NULL);
	assumpB.SetImplementsLinkAPI(true);

	NxsCharactersBlock charsB(NULL, NULL);
	charsB.SetCreateImpliedBlock(true);
	charsB.SetImplementsLinkAPI(true);
	charsB.SetSupportMixedDatatype(true);
	charsB.SetConvertAugmentedToMixed(true);
	
	NxsDataBlock dataB(NULL, NULL);
	dataB.SetCreateImpliedBlock(true);
	dataB.SetImplementsLinkAPI(true);
	dataB.SetSupportMixedDatatype(true);
	dataB.SetConvertAugmentedToMixed(true);
	
	NxsTaxaBlock taxaB;
	taxaB.SetImplementsLinkAPI(false);
	
	NxsTreesBlock treesB(NULL);
	treesB.SetCreateImpliedBlock(true);
	treesB.SetImplementsLinkAPI(true);
	treesB.SetProcessAllTreesDuringParse(true);
	treesB.SetAllowImplicitNames(true);

	std::string emptyString;
	NxsStoreTokensBlockReader  storerB(emptyString, true);
	storerB.SetImplementsLinkAPI(false);


	if (!treesBlockOnly) {
		factory.AddPrototype(&assumpB, "ASSUMPTIONS");
		factory.AddPrototype(&assumpB, "SETS");
		factory.AddPrototype(&assumpB, "CODONS");
		factory.AddPrototype(&charsB, "CHARACTERS");
		factory.AddPrototype(&dataB, "DATA");
		factory.AddPrototype(&taxaB);
	}
	factory.AddPrototype(&treesB);
	nexusReader.Add(this);

	try {
		ioObject.message << "Executing " << filepath;
		ioObject.printMessage(BullIO::STATUS_MSG_LEVEL);
		readerStack.push(&nexusReader);
		nexusReader.Execute(token);
	}
	catch(...) {
		nexusReader.RemoveFactory(&factory);
		nexusReader.Detach(this);
		readerStack.pop();
		BlockReaderList blocks = nexusReader.GetUsedBlocksInOrder();
		for (BlockReaderList::iterator b = blocks.begin(); b != blocks.end(); ++b) {
			if (*b != this)
				delete *b;
		}
		throw;
	}

	readerStack.pop();
	nexusReader.RemoveFactory(&factory);
	nexusReader.Detach(this);
	return nexusReader.GetUsedBlocksInOrder();
}

bool BullShell::readFile(const std::string &filepath, const bool treesBlockOnly) 
{
	try {
		BlockReaderList blocks = getBlocksFromFile(filepath, treesBlockOnly);

		this->identifyAndProcessUnhandledBlocks(blocks); /* deal with all of the non-BULL blocks */
		
		for (BlockReaderList::iterator b = blocks.begin(); b != blocks.end(); ++b) {
			if (*b != this)
				delete *b;
		}
		return true;
	}
	catch (const NxsException &x) {
		NexusError(x.msg, x.pos, x.line, x.col);
		return false;
	}
	return true;
}

void BullShell::handleOneLineCommand(std::string cmd)
{
	std::string fullBlock;
	fullBlock.append("; ");
	fullBlock.append(cmd);
	fullBlock.append(" ; end; ");
	istringstream cmdin(fullBlock.c_str());

	NexusToken token(cmdin);
	try {
		Read( token );
	}
	catch(const NxsException & x ) {
		if (errormsg.empty())
			NexusError(x.msg, x.pos, x.line, x.col );
		else	
			NexusError(errormsg, x.pos, x.line, x.col );
	}
}

void BullShell::Run(const char *infn)
{
	quit_now = false;
	try {
		std::string cmd;
		if (infn != NULL) {
			cmd.assign("Execute ");
			std::string sfn(infn);
			cmd.append(NxsString::GetEscaped(sfn));
			handleOneLineCommand(cmd);
		}
		while ( !quit_now ) {
			cerr << endl;
			cerr << "BullShell> ";
			cmd.clear();
			cin >> cmd;
			if (cin.eof())
				return;
			handleOneLineCommand(cmd);
		}
	}
	catch(bad_alloc &ba) {
		cout << "There is not enough memory to execute that command, BullShell will be crashing now."<< std::endl;
		throw;
	}	
}



/**
 * @method Read [void:protected]
 * @param token [NexusToken&] the token used to read from in
 * @throws XBull
 *
 * This function provides the ability to read everything following
 * the block name (which is read by the Nexus object) to the end or
 * endblock statement. Characters are read from the input stream
 * in. Overrides the pure virtual function in the base class.
 */
void BullShell::Read( NexusToken& token )
{	
	processInterveningPublicNexusBlocks();
	isEmpty = false;
	isUserSupplied = true;
	DemandEndSemicolon(token, "Begin BULL");

	for(; ; ) {
		token.GetNextToken();
		NxsBlock::NxsCommandResult res = HandleBasicBlockCommands(token);
		if (res == NxsBlock::NxsCommandResult(STOP_PARSING_BLOCK))
			return;
		if (res != NxsBlock::NxsCommandResult(HANDLED_COMMAND)) {
			if (token.Equals("Quit"))
				{
				quit_now = true;
				ioObject.message << "\nBull says goodbye\n"; 
				ioObject.printMessage(BullIO::INFO_MSG_LEVEL); 
				break;
				}
			else if (token.Abbreviation("CODLIKESTartval") || token.Abbreviation("CLIKESTartval") || token.Abbreviation("CLSTartval"))
				HandleCodLikeStartVal(token );
			else if ( token.Abbreviation("EXecute")  ) 
				HandleExecute(token);
			else if ( token.Abbreviation("Gettrees")) 
				HandleGetTrees( token );		
			else if ( token.Abbreviation("LOg") ) 
				HandleLog( token ); 
			else if ( token.Abbreviation("SImulate") ) 
				HandleSimulate( token );	
			else
				SkipCommand(token);
		}
	}	
}


/**
 * @method NexusError [virtual void:public]
 * @param msg [std::string&] the error message
 * @param pos [streampos] the point in the NEXUS file where the error occurred
 * @param line [long] the line in the NEXUS file where the error occurred
 * @param col [long] the column in the NEXUS file where the error occurred
 *
 * Called when an error is encountered in a NEXUS file. Allows program to
 * give user details of the error as well as the precise location of the
 * error. Virtual function that overrides the pure virtual function in the
 * base class Nexus.
 */
void BullShell::NexusError(
	const string& msg, 
	file_pos  pos, 
	long line, 
	long col)
{
	ioObject.message << msg << '\n';
	if (line >=0)
		ioObject.message << "at line " << line << ", column (approximately) " << col << " (and file position "<< int(pos) << ")\n";
	ioObject.printMessage(BullIO::ERROR_MSG_LEVEL);
}


/**
 * @method Report [virtual void:public]
 * @param out [std::ostream&] the output stream to which to write the report
 *
 * This function outputs a brief report of the contents of this BullShell block.
 * Overrides the pure virtual function in the base class.
 */
void BullShell::Report( std::ostream& /* out */ )
{}


/**
 * @method FindTreeFromName [Tree *:protected]
 * @param s [std::string] the name of the requested tree
 *
 * searches through the vector of Tree * for the requested tree
 * not case sensitive because Tree names are capitalized
 * NOTE it doesn't check to see if there is more than one Tree in the
 * list with the correct name, it just returns the first.
*/
Tree * BullShell::FindTreeFromName(std::string s) {
	return kernel.FindTreeFromName(s);
}

bool bull::FileExists( const char* fn )
{	
	ifstream testst;
	testst.open(fn);
	if (testst.good()) {
		testst.close();
		return true;
	}
	testst.close();
	return false;
}
