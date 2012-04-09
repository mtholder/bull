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
 

#include <algorithm>
#include "settings.hpp"

using namespace bull;


LikePartSettings::LikePartSettings(
  const PartSettings &ps,
  const std::vector<Model *> & modVec) 
  :sharedPartSettings(ps)
{
	modelMixingSetting=par(CUR);
	for (std::vector<Model *>::const_iterator i = modVec.begin(); i != modVec.end(); ++i)
		modSettings.push_back(ModSettings(*i));
}

LikePartSettings::~LikePartSettings() {
}


void LikePartSettings::SetModel(Model *m, std::size_t n) {
	return modSettings.at(n).SetModel(m);
}


void SimSettings::AddPartSetting(PartSettings ps, const std::vector<Model *> & modVec) {
	partSettingVec.push_back(new LikePartSettings(ps, modVec));
}

void SimSettings::SetOutputDatatype(std::size_t otIndex,int o) 
{	
	if (this->outputType.size() > otIndex) {
		this->outputType[otIndex] = o;
		return;
	}
	assert(this->outputType.size() == otIndex); //must be called inorder for the push_back to work
	this->outputType.push_back(o);
}


/* reset all nonpersistent fields*/
void SimSettings::Reset() {
	appToOut = false;
	concats = 1;
	nReps = 1;
	nSimsAtOnce = 1;
	numOutputs = 1;
	outputType.clear();

	std::vector<LikePartSettings *>::iterator psIt = partSettingVec.begin();
	for (; psIt != partSettingVec.end(); ++psIt)
		delete *psIt;
	partSettingVec.clear();
}

/// assumes that `this` is empty (Reset has been called or fresh instance)
void SimSettings::copyData(const SimSettings & other) {
	this->appToOut = other.appToOut;
	this->concats = other.concats;
	this->nReps = other.nReps;
	this->nSimsAtOnce = other.nSimsAtOnce;
	this->outputType = other.outputType;

	assert(this->partSettingVec.empty());

	std::vector<LikePartSettings *>::const_iterator psIt = other.partSettingVec.begin();
	for (; psIt != other.partSettingVec.end(); ++psIt)
		this->partSettingVec.push_back(new LikePartSettings(**psIt));
	
}

unsigned PartSettings::GetRawDataSize() const {
	unsigned s = 0;
	std::vector<segment> ::const_iterator seg = this->segmentsInPartition.begin();
	for (; seg != this->segmentsInPartition.end(); ++seg)
		s += (unsigned)(seg->second - seg->first); 
	return s;
}

ModSettings::ModSettings(Model *m)
{	useCurrentBrLen=true;
	stationaryModel=true;
	startingBrLen=bull::DEFAULT_START_EDGE_LEN;	
	mod=m;
	blenModSetting=blenSetting=par(DEF)|par(MIN);
	hasBrLenMod=false;
	numBrLenMod=0;
}



