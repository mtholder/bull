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
 


#include "string_extensions.hpp"


std::string	   &AppendDouble(std::string& s, const double d)
{
	char tmp[81];
	sprintf( tmp, "%.6f", d );
	int tmplen = strlen(tmp);
	for (; ;) {
		if( tmplen < 3 || tmp[tmplen-1] != '0' || tmp[tmplen-2] == '.' )
			break;
		tmp[tmplen-1] = '\0';
		tmplen--;
	}
	s.append(tmp);
	return s;
}


std::string GetTreeToken(std::string fullstr,std::string::iterator &c)
{	//Expecting either a name or tree code character (; ,:) 
	//Names: can be in single quotes (which are stripped off).	
	//Preceeding and terminal whitespace is removed (except if it is in a single quoted name)
	
	std::string tk;
	const char ch = GetNextGraphicalNonCommentedChar(fullstr,c);
	AppendChar(tk, ch);
	if(*c == '\'')	
			{tk="";
			while (*c!='\'' && c!=fullstr.end())
				tk+=*c++;
			if(c!=fullstr.end())	c++;
			return tk;
			}
	if(c == fullstr.end())	return tk;
	c++;
	if(tk[0] == '(' || tk[0] == ')' || tk[0] == ':' || tk[0] == ';' || tk[0] == ','|| *c == '[')
		return tk;
	while (*c!='(' && *c!=')' && *c!=':' && *c!=';' && *c!=','&& *c!='\'' && isgraph(*c))
		{if(*c == '[')
			GetNextNonCommentedChar(fullstr,c);
		else
			{tk+=*c;
			c++;
			}
		}
	return tk;	
}

char GetNextGraphicalNonCommentedChar(std::string fullstr,std::string::iterator &c)
{	while (c!=fullstr.end() &&!isgraph(*c) )	c++;
	if(c == fullstr.end())	return 0;
	if(*c!='[') return *c;
	while (c!=fullstr.end() && *c!=']') c++;
	if(c == fullstr.end())	return 0;
	c++;
	return GetNextGraphicalNonCommentedChar(fullstr,c); 
}

char GetNextNonCommentedChar(std::string fullstr,std::string::iterator &c)
{	if(*c!='[') return *c;
	while (c!=fullstr.end() && *c!=']') c++;
	if(c == fullstr.end())	return 0;
	c++;
	return GetNextNonCommentedChar(fullstr,c);	
}

std::string& BlanksToUnderscores( std::string& s )
{
   int len = s.length();
   for ( int k = 0; k < len; k++ ) {
	  if( s[k] == ' ' )
		 s[k] = '_';
   }
   return s;
}

std::string& UnderscoresToBlanks( std::string& s )
{
   int len = s.length();
   for ( int k = 0; k < len; k++ ) {
	  if( s[k] == '_' )
		 s[k] = ' ';
   }
   return s;
}

//breaks the input string into words as read by >> and returns the number of words
int BreakString(std::vector<std::string> &words,std::string phrase)
{	int nwords=0;
	int l=phrase.length();
	//Count words in phrase;
	bool grap=false;
	int startofThisWord=-1;
	for (int i=0; i<l; i++)	
		{if(isgraph(phrase[i]))
			{if(!grap) 
				{nwords++;
				grap=true;
				startofThisWord=i;
				}
			}
		else
			{
			if(grap)
				{grap=false;
				words.push_back(std::string(phrase,startofThisWord,i-startofThisWord));
				startofThisWord=-1;
				}
			}
		}
	if(startofThisWord>-1)
		words.push_back(std::string(phrase,startofThisWord,l-startofThisWord));
	return nwords;
	
}


