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
 

#ifndef LIKESETTINGS
#define LIKESETTINGS
#include <map>

#include <cstddef>
#include <utility>
#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 

#include "basic_bull.hpp"
#include "ssrf_codon_sub_mod.hpp"
#include "mth_exception.hpp"
#include "ncl/nxsdefs.h"
#include "part_mod_index.hpp"
namespace bull {

class LikeSettingException: public MTHException
{
	public:
		LikeSettingException(const char *p)
			:MTHException(p)
			{}
};

class ModSettings
{
	public :
		ModSettings(Model *);
		
		const Model *GetModel() const {
			return this->mod;
		}
		Model *GetModel() {
			return this->mod;
		}
		void SetModel(Model *m) {
			this->mod = m;
		}
		void SetBrLenOwner(PartModIndex pm) {
			this->ownerOfBlen = pm;
		}
		PartModIndex GetBrLenOwner() const {
			return this->ownerOfBlen;
		}
		PartModIndex GetBrLenModOwner() const {
			return this->ownerOfBlenMod;
		}
		int GetBLenSetting() const {
			return this->blenSetting;
		}
		int GetBLenModSetting() const {
			return this->blenModSetting;
		}
		bool HasBrLenMod() const {
			return this->hasBrLenMod;
		}
		unsigned GetNumBrLenMod() const {
			return this->numBrLenMod;
		}
	protected:
		bool stationaryModel;
		bool useCurrentBrLen;
		double startingBrLen;
		Model *mod;

		// change ownerOfBlenMod to vector to allow relaxed clock starting 
		//	points (multiple modifiers with one model)
		PartModIndex ownerOfBlen;
		PartModIndex ownerOfBlenMod;
		int blenSetting;
		int blenModSetting;
		unsigned numBrLenMod;
		bool hasBrLenMod;
};

class PartSettings
{
	public :
		typedef std::pair<unsigned, unsigned> segment;

		
		PartSettings(unsigned start,unsigned end);

		unsigned GetRawDataSize() const;

		std::vector<segment> segmentsInPartition;

};

inline PartSettings::PartSettings(unsigned s, unsigned e)
{
	assert(s <= e);
	segmentsInPartition.push_back(segment(s,e)); 
}

class LikePartSettings
{
	public:
		LikePartSettings(const PartSettings &, const std::vector<Model *> &);
		~LikePartSettings();

		bool AllStationaryModels() const;
		const Model *GetModel(std::size_t i) const {
			return modSettings[i].GetModel();
		}
		Model *GetModel(std::size_t i) {
			return modSettings[i].GetModel();
		}
		int GetModelMixingParamSettings() const {
			return this->modelMixingSetting;
		}
		double GetStartingFreqnOfThisModel(std::size_t j) const;
		PartModIndex GetBrLenOwner(std::size_t j) const {
			return this->modSettings.at(j).GetBrLenOwner();
		}
		PartModIndex GetBrLenModOwner(std::size_t j) const {
			return this->modSettings.at(j).GetBrLenModOwner();
		}
		int GetBLenSetting(std::size_t i) const {
			return this->modSettings.at(i).GetBLenSetting();
		}
		int GetBLenModSetting(std::size_t i) const {
			return this->modSettings.at(i).GetBLenModSetting();
		}
		bool HasBrLenMod(std::size_t j) const {
			return this->modSettings.at(j).HasBrLenMod();
		}
		int GetNumBrLenMod(std::size_t m) const {
			return this->modSettings.at(m).GetNumBrLenMod();
		}

		void SetUseCurrentBrLen(bool);
		void SetModel(Model *m, std::size_t i);
		void SetBrLenOwner(std::size_t j, PartModIndex pmi) {
			this->modSettings.at(j).SetBrLenOwner(pmi);
		}
		
		const PartSettings & GetPartSettings() const {
			return this->sharedPartSettings;
		}
		
	protected:
		std::vector<ModSettings> modSettings;
		int modelMixingSetting;
		PartSettings sharedPartSettings;
		
};

class SimSettings {
	public:
		SimSettings() {
			Reset();
		}
		~SimSettings() {
			Reset();
		}

		bool getAppendToOut() const {
			return this->appToOut;
		}
		unsigned GetNReps() const {
			return nReps;
		}
		unsigned GetNDataPartitions() const {
			return this->partSettingVec.size();
		}
		unsigned GetPartitionLength(std::size_t i) const {
			assert(this->partSettingVec.size() > i);
			return this->partSettingVec[i]->GetPartSettings().GetRawDataSize();
		}
		Model * GetModel(std::size_t i,std::size_t j) {
			assert(this->partSettingVec.size() > i); 
			return this->partSettingVec[i]->GetModel(j);
		}
		unsigned getNumConcats() const {
			return this->concats;
		}
		unsigned getNumSimsAtOnce() const {
			return this->nSimsAtOnce;
		}
		int GetOutputDatatype(std::size_t otIndex) const {
			assert(otIndex < this->outputType.size());
			return this->outputType[otIndex];
		}
		
		void AddPartSetting(PartSettings ps, const std::vector<Model *> &);

		void setAppendToOut(bool a) {
			this->appToOut = a;
		}
		void SetOutputDatatype(std::size_t otIndex, int o);
		void setNumConcats(unsigned c) {
			this->concats = c;
		}
		void setNumOutputs(unsigned n) {
			this->numOutputs = n;
		}
		void setNumSimsAtOnce(unsigned n) {
			this->nSimsAtOnce = n;
		}
		void setNumReps(unsigned n) {
			this->nReps = n;
		}
		void Reset();
		
		SimSettings(const SimSettings &other) {
			*this = other;
		}
		SimSettings & operator=(const SimSettings & other) {
			this->Reset();
			this->copyData(other);
			return *this;
		}

	protected:
		void copyData(const SimSettings & other);
		bool appToOut;
		unsigned concats;
		unsigned nReps;
		unsigned nSimsAtOnce;
		unsigned numOutputs;
		std::vector<LikePartSettings *> partSettingVec;
		std::vector<int> outputType;
	private:
};


} // namespace bull

#endif

