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
 

#if !defined(BULLCMDLINE_H)
#define BULLCMDLINE_H

#include <ostream>
#include <stack>


#include "ncl/nxsblock.h"
#include "ncl/nxsreader.h"
class NxsAssumptionsBlock;
class NxsCharactersBlockAPI;
class NxsTaxaBlockAPI;
class NxsTreesBlockAPI;

#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 
#if defined HAVE_UNISTD_H && HAVE_UNISTD_H
#	include <unistd.h>
#elif defined HAVE_SYS_TYPES_H && HAVE_SYS_TYPES_H
#	include <sys/types.h>
#endif


#include "bull_kernel.hpp"

#include "basic_bull.hpp"
#include "bull_io.hpp"
#include "bull_listener.hpp"
#include "bull_parsers.hpp"
#include "char_encoding.hpp"
#include "encoded_chars.hpp"
#include "util.hpp"
#include "tree.hpp"
#include "ssrf_codon_sub_mod.hpp"
#include "settings.hpp"

namespace bull {
class NoSuchTree : public MTHException 
{
	public :
		NoSuchTree() : MTHException() {}
		NoSuchTree(const char *c) :MTHException(c) { }
};

class NoModel : public MTHException { public :
	NoModel() : MTHException() {}
	NoModel(const char *c) :MTHException(c) { }
	};

class BullKernel;
class BullShell:  public NxsBlock, public BullListener
{
	public:
		BullShell(BullKernel & kernel);

		void OutputComment( std::string s );
		void ExitingBlock( std::string blockName );
		void EnteringBlock( std::string blockName );
		void NexusError(const std::string& msg, file_pos pos, long line, long col );
		void Report(std::ostream& out );
		void Run(const char *);

		void stateChanged(Event, void * , void *) {
		}
	private:
		void identifyAndProcessUnhandledBlocks(std::list<NxsBlock*> &);
		BlockReaderList getBlocksFromFile(const std::string &filepath, const bool treesBlockOnly);
		void handleOneLineCommand(std::string cmd);
		void processInterveningPublicNexusBlocks();
		void processNexusBlocks(std::list<NxsBlock*> &);
		
		void processAssumptionsBlock(NxsAssumptionsBlock *);
		void processCharactersBlock(NxsCharactersBlockAPI *);
		void processTaxaBlock(NxsTaxaBlockAPI *);
		void processTreesBlock(NxsTreesBlockAPI *);
		bool readFile(const std::string &filepath, const bool treesBlockOnly = false);
		
		Tree *FindTreeFromName(std::string);

		EncodingType InterpretOutputType(NexusToken &token) const;
		void RequireSemicolon(NexusToken & token, const char *cmd) const;
		void ThrowNoSemicolon(NexusToken & token, const char *cmd) const;

		void HandleEndblock( NexusToken& token );
		void HandleExecute( NexusToken& token );
		void HandleExecuteBullBlocksInFile( NexusToken& tok );
		void HandleLog( NexusToken& token );
		void HandleGetTrees(NexusToken& token );
		void HandleCodLikeStartVal(NexusToken& token);
		void HandleSimulate(NexusToken& token );

		
		GetTreesOpts ParseGetTrees(NexusToken& ttoken) const;
		void ParseCodLikeStartValCommand(CodLikeStartOpts &, NexusToken& ttoken) const;
		SimulateOpts ParseSimulateCommand(NexusToken& token) const;
		
		void ExecGetTrees(const GetTreesOpts & gto);
		void ExecCodLikeStartValCommand(const CodLikeStartOpts &);
		void ExecSimulateCommand(const SimulateOpts &so)
		{
			kernel.simulate(so);
		}

		void Read( NexusToken& token );
		void Reset()
			{
			kernel.Reset();
			ioObject.Reset();
			}


		mutable BullIO ioObject;
		bool quit_now;
		BullKernel & kernel;
		std::stack<NxsReader *> readerStack;
};

} //namespace bull 


#endif


