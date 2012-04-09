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
#include <iostream>
#include "bull_io.hpp"
using namespace std;
using namespace bull;

void BullIO::stopLogging()
{
	if (logf_open) {
		logf.close();
		logf_open = false;
		message = "Logging to file stopped.\n";
		printMessage(STATUS_MSG_LEVEL);
	}
}

void BullIO::startLogging(const std::string logfname, bool append)
{
	if (logf_open)
		stopLogging();
	const bool exists = FileExists( logfname.c_str() );
	
	if (append) {
		logf_open = true;
		logf.open( logfname.c_str(), ios::out | ios::app );
		if (exists) {
			message << "\nAppending to log file " << logfname;
			printMessage(STATUS_MSG_LEVEL);
		}
	}
	else {
		logf_open = true;
		logf.open( logfname.c_str() );
		if (exists) {
			message << "\nReplacing log file " << logfname;
			printMessage(STATUS_MSG_LEVEL);
		}
	}
	if (!exists) {
		message << "\nCreating new log file " << logfname;
		printMessage(STATUS_MSG_LEVEL);
	}
}

BullIO::BullIO()
{
	Reset();
}

BullIO::~BullIO() 
{
	Reset();
}

void BullIO::Reset() 
{
	min_visible_level = STATUS_MSG_LEVEL;
	if (logf_open)
		logf.close();
	inf_open = false;
	logf_open = false;
}

/**
 * Asks user if "something" is ok, where "something" is expressed
 * in the title and message displayed.	This is a virtual function
 * so it can be overridden in a derived class to use a different
 * (perhaps graphical) means for displaying the message.
 * Note: mb_message should terminate with a quesiton mark; none
 * will be provided by this function.
 */
bool BullIO::userSaysOk(const std::string & mb_message, const std::string & mb_title, bool defaultVal) const 
{
	const unsigned MAX_RESPONSE_LEN = 255;
	char response[MAX_RESPONSE_LEN];
	
	while (true) {
		cout << '\n' << mb_title << "\n  " << mb_message << " (y/n) ";	
		cin.getline(response, MAX_RESPONSE_LEN);

		if (!cin.good()) {
			if (cin.eof())
				return defaultVal;
			cin.clear();
		}
		else if (response[1] == '\0') {
			char c = response[0];
			if (c == 'y' || c == 'Y')
				return true;
			if (c == 'n' || c == 'N')
				return false;
		}
		cout << "\nMust answer by typing either y or n and then pressing the Enter key\n\n";
	}
}

static const char * MSG_PREFIXES[6] = {"DEBUG: ", "", "", "Warning: ", "Error: ", "Fatal Error: "};
void BullIO::print(
  const MessageLevel ml,
  const std::string & msg,
  bool linefeed, 
  bool flusht) const 
{
  	int iml  = (int) ml;
  	if (iml < (int) min_visible_level)
  		return;
  	
  	if (iml < 0 || iml > 5)
  		iml = 2;
  	const char *pref = MSG_PREFIXES[iml];
  	cout << pref << msg; 
	if ( linefeed ) {
		cout << '\n';
		if (flusht)
			cout << flush;
	}	
	if (logf_open) {
		logf << msg;
		if ( linefeed ) {
			logf << '\n';
			if (flusht)
				logf << flush;
		}
	}
}

