// This file is part of BULL, a program for phylogenetic simulations
// most of the code was written by Mark T. Holder.

//It is distributed in the hope that it will be useful,
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

#if ! defined(PART_MOD_INDEX_HPP)
#define PART_MOD_INDEX_HPP

#include <cstddef>

namespace bull {

///Class that serves holds an index to a partition subset and an index to model 
///	within the subset.

class PartModIndex
{
	public:
		std::size_t partNum;
		std::size_t modNum;

		PartModIndex()
			:partNum(0),
			modNum(0)
			{}
		PartModIndex(std::size_t p, std::size_t m)
			:partNum(p),
			modNum(m)
			{}

		bool operator<(const PartModIndex s) const {
			if (this->partNum == s.partNum)
				return this->modNum < s.modNum;
			return this->partNum < s.partNum;
		}
		bool operator!=(const PartModIndex s) const {
			return (this-> partNum != s.partNum) || (this->modNum != s.modNum);
		}
		bool operator == (const PartModIndex s) const {
			return (this->partNum == s.partNum) && (this->modNum == s.modNum);
		}

		void Set(std::size_t p,std::size_t m) {
			this->partNum = p;
			this->modNum = m;
		}
};

} // namespace bull

#endif
