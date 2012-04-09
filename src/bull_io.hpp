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
#if !defined(BULL_IO_HPP)
#define BULL_IO_HPP

#include <fstream>
#include <string>

#include "ncl/nxsstring.h"

namespace bull {

bool FileExists( const char* fn );
class BullIO
{
	public:
		enum MessageLevel {
			DEBUG_MSG_LEVEL,
			STATUS_MSG_LEVEL,
			INFO_MSG_LEVEL,
			WARN_MSG_LEVEL,
			ERROR_MSG_LEVEL,
			FATAL_MSG_LEVEL
			};

		BullIO();
		virtual ~BullIO();

		void Reset();

		virtual void	print(const MessageLevel, const std::string &msg, bool linefeed = true ,bool flusht = true) const;
		void 			printMessage(MessageLevel ml, bool linefeed = true , bool flusht = true) {
			this->print(ml, this->message, linefeed, flusht);
			this->message.clear();
		}
		virtual bool userSaysOk(const std::string & mb_message, const std::string & mb_title, bool defaultVal) const;

		NxsString message;
		
		MessageLevel GetFilterLevel() const
			{
			return min_visible_level;
			}

		void SetFilterLevel(MessageLevel ml)
			{
			min_visible_level = ml;
			}
		void stopLogging();
		void startLogging(const std::string logfname, bool append);

	protected:
		bool inf_open;
		bool logf_open;
		mutable std::ofstream logf;
		MessageLevel min_visible_level;
};

} //namespace bull 

#endif
