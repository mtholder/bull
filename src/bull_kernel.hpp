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
#if !defined(BULL_KERNEL_HPP)
#define BULL_KERNEL_HPP

#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <set>

#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 
#include "basic_bull.hpp"
#include "settings.hpp"
#include "bull_listener.hpp"
#include "ssrf_codon_sub_mod.hpp"
namespace bull {

extern long gRngSeed;

std::string EncodingTypeToString(EncodingType et);

class CodLikeStartOpts;
class SetOfLikeAttr;
class SetOfTreeSimAttr;
class SimulateOpts;
class Tree;

enum BullAttribute {
	kChars         = 0x01,
	kModel         = 0x02,
	kCharMod       = 0x03,
	kTaxa          = 0x04,
	kCharTax       = 0x05,
	kModTax        = 0x06,
	kCharModTax    = 0x07,
	kTrees         = 0x08,
	kCharTre       = 0x09,
	kModelTre      = 0x0A,
	kCharModTre    = 0x0B,
	kTaxTre        = 0x0C,
	kCharTaxTre    = 0x0D,
	kModTaxTre     = 0x0E,
	kCharModTaxTre = 0x0F
	};

/// The update mode is modelled after the mode of PAUP's gettrees.  The enum
/// uses bit fields to refer to elements that are:
///		unique to the new data
///		unique to the old data
///		found in both new and old. 
/// Using this enum one can specify what the behavior with respect to new data.
enum UpdateMode {				
	kClearUpdateMode     = 0x00, /// stored = {}
	kOnlyNewUpdateMode   = 0x01, /// stored = new - {stored&new}
	kIntersectUpdateMode = 0x02, /// stored = {stored&new}
	kReplaceUpdateMode   = 0x03, /// stored = new
	kOnlyOldUpdateMode   = 0x04, /// stored = stored - {stored&new}
	kXorUpdateMode       = 0x05, /// stored = {stored^new}
	kIgnoreUpdateMode    = 0x06, /// stored = stored
	kAppendUpdateMode    = 0x07  /// stored = stored|new
	};

class BullKernel
{
	public:
		BullKernel();
		~BullKernel();
	
		void addListener(BullListener *listener) {
			this->listeners.insert(listener);
		}
		void removeListener(BullListener *listener) {
			this->listeners.erase(listener);
		}
		
		unsigned getNumTaxa() const {
			return taxaLabels.size();
		}

		unsigned getNumTrees() const {
			return treeList.size();
		}
		
		Tree *getTree(unsigned i) const {
			if (i >= treeList.size())
				return NULL;
			std::list<Tree *>::const_iterator trIt =  treeList.begin();
			for (; trIt != treeList.end(); ++trIt, --i) {
				if (i == 0)
					return *trIt;
			}
			return NULL;
		}
		
		bool hasModel() const;
		
		const std::vector<std::string> & updateTaxa(const std::vector<std::string> & taxa, const UpdateMode mode = kAppendUpdateMode);
		std::vector<Tree *> updateTrees(std::vector<Tree *> trees, const UpdateMode mode = kAppendUpdateMode);
		void updateModels(const SimSettings &, const UpdateMode mode = kAppendUpdateMode);

		void simulate(const SimulateOpts &);

		Tree * FindTreeFromName(std::string s);

		void Reset();

		static void WriteSimulations(unsigned num, Tree *t, std::string outpfn,
									 unsigned outputEncodingType, unsigned totRepNum, 
									 long seedAtBeg, unsigned simCondNum,  const SimSettings &);
		void setModelConstRef(const SSRFCodonSubModSet &ms) {
			modelSet = ms;
		}
		const SSRFCodonSubModSet & getModelConstRef() const {
			return modelSet;
		}
		
	private:
		static SetOfLikeAttr *CreateSetOfSimAttrForANode(SimSettings & sSettings);
		static SetOfTreeSimAttr *CreateSetOfSimAttrForATree(SimSettings & sSettings);
		static void DestroySimulationAttributes(Tree &t);
		static void PrepareTreeForSimulation(Tree &, SimSettings & sSettings);
		
		void deleteModels(bool leaveOpen);
		void deleteTrees(bool leaveOpen);


		void openChars() {
			if (charsOpen)
				return;
			charsOpen = true; 
			notifyListeners(BullListener::kCharsOpened);
		}
		void closeChars() {
			assert (charsOpen);
			charsOpen = false;
			notifyListeners(BullListener::kCharsClosed);
		}
		void openModels() {
			if (modelsOpen)
				return;
			modelsOpen = true;
			notifyListeners(BullListener::kModelsOpened);
		}
		void closeModels() {
			assert (modelsOpen);
			modelsOpen = false;
			notifyListeners(BullListener::kModelsClosed);
		}
		void openTrees() {
			if (treesOpen)
				return;
			treesOpen = true;
			notifyListeners(BullListener::kTreesOpened);
		}
		void closeTrees() {
			assert(treesOpen);
			treesOpen = false;
			notifyListeners(BullListener::kTreesClosed);
		}
		void openTaxa() {
			if (taxaOpen)
				return;
			taxaOpen = true;
			notifyListeners(BullListener::kTaxaOpened);
		}
		void closeTaxa() {
			assert(taxaOpen);
			taxaOpen = false;
			notifyListeners(BullListener::kTaxaClosed);
		}
		void notifyListeners(BullListener::Event, void *blob = 0L);


		Criterion criterion;
		int datatype;
		AnalysisMode analysisMode;

		std::vector<std::string> taxaLabels;
		std::map<std::string, unsigned> taxLabelToInd;
		
		std::list<Tree *> treeList; /// trees in memory (use deleteTrees() to free)
		
		int numSimulationConditionsDone;

		bool workingOnBranches;
		SSRFCodonSubModSet modelSet;

		
		bool charsOpen;
		bool modelsOpen;
		bool taxaOpen;
		bool treesOpen;
		std::set<BullListener*> listeners;
		
		friend class BullAttributeManager;

	public:
		bool logEachStep;		
};


/// friend of the kernel that creates exception safe open and close statements
/// 
class BullAttributeManager
{
	public:
		~BullAttributeManager() {
			if (att & BullAttribute(kChars))
				kern.closeChars();
			if (att & BullAttribute(kModel))
				kern.closeModels();
			if (att & BullAttribute(kTaxa))
				kern.closeTaxa();
			if (att & BullAttribute(kTrees))
				kern.closeTrees();
		}
	private:
		friend class BullKernel;
		BullAttributeManager(BullAttribute a, BullKernel & kernel)
			:att(a),
			kern(kernel) {
			if (att & BullAttribute(kChars))
				kern.openChars();
			if (att & BullAttribute(kModel))
				kern.openModels();
			if (att & BullAttribute(kTaxa))
				kern.openTaxa();
			if (att & BullAttribute(kTrees))
				kern.openTrees();
		}
			
		const BullAttribute att;
		BullKernel & kern;
};

} //namespace bull 

#endif
