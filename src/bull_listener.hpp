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
#if !defined(BULL_LISTENER_HPP)
#define BULL_LISTENER_HPP

namespace bull {

class BullListener 
{
	public:
		enum Event {
			kInitialized,
			kDying,
			kDead,
			kCharsOpened,
			kCharsDeleted,
			kCharsAdded,
			kCharsClosed,
			kModelsOpened,
			kModelsDeleted,
			kModelsAdded,
			kModelsClosed,
			kTaxaOpened,
			kTaxaDeleted,
			kTaxaAdded,
			kTaxaClosed,
			kTreesOpened,
			kTreesDeleted,
			kTreesAdded,
			kTreesClosed,
			kSimStarted,
			kSimAvailable,
			kSimFinished
		};
		virtual ~BullListener() {}
		virtual void stateChanged(Event, void *target, void *blob) = 0;
};


} //namespace bull 

#endif
