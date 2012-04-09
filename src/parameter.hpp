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
 

#ifndef PARAMETER
#define PARAMETER


#if defined HAVE_CONFIG_H && HAVE_CONFIG_H
#	include <config.h>
#endif 

#include "basic_bull.hpp"
#include "mth_exception.hpp"

namespace bull {


class ParamExcep: public MTHException
{
	public: 
		ParamExcep(const char *p)
			:MTHException(p) 
			{}
		ParamExcep(std::string n, const char *p)
			:MTHException(n,p) 
			{}
};


class ParamOutOfRangeExcep: public MTHException 
{
	public:
		ParamOutOfRangeExcep(const char *p)
			:MTHException(p) 
			{}
		ParamOutOfRangeExcep(std::string n, const char *p)
			:MTHException(n,p) 
			{}
		ParamOutOfRangeExcep() 
			:MTHException() 
			{}
};

//DEF, CUR , RAN , APPRO and APPRO2 are starting conditions, make sure that ALLSTART is their union if adding another starting poing
enum par {
	MIN=1, 
	MAX=2, 
	DEF=4, 
	CUR=8, 
	RAN=16, 
	APPRO = 32, 
	FIX = 64, 
	UBOUN = 128, 
	LBOUN = 256, 
	CONSTR = 512,
	APPRO2 = 1024, 
	ALLSTART = 1084, 
	CONSTR2 =2048
	};

class Parameter 
{
	protected :
		double defaultv;
		int setting;
	public:
		double val;
		Parameter();
		Parameter(double v,int st,double def);
		Parameter(double v,int st);
		Parameter(double v);
		virtual ~Parameter(){
		}
		virtual bool HasMin() const {
			return false;
		}
		virtual void SetMin(double ) {
			throw ParamExcep("Calling inappropriate base function of Parameter");
		}
		virtual bool HasMax() const {
			return false;
		}
		virtual void SetMax(double) {
			throw ParamExcep("Calling inappropriate base function of Parameter");
		}
		int StartWithDefault();
		virtual void SetToDefault();
		virtual void SetDefault(double m);
		virtual void SetStartWithDefault(bool i);
		int StartWithCurrent()	;
		virtual void SetCurrent(double m)	{val=m; }
		virtual void SetStartWithCurrent(bool i);	
		int StartWithRandom();
		virtual void SetStartWithRandom(bool i) ;
		int StartWithApproximation()	;
		virtual void SetStartWithApproximation(bool i);
		int Fixed();
		virtual void SetFixed(double m) ;
		void SetFixed(bool i);
		virtual bool HasUpperBound() const {
			return false;
		}
		virtual void SetUpperBound(double ) {
			throw ParamExcep("Calling inappropriate base function of Parameter");
		}
		virtual void SetUpperBound(bool ) {
			throw ParamExcep("Calling inappropriate base function of Parameter");
		}
		virtual bool HasLowerBound() {
			return false;
		}
		virtual void SetLowerBound(double m)	{throw ParamExcep("Calling inappropriate base function of Parameter"); m=1; }
		virtual void SetLowerBound(bool i)		{throw ParamExcep("Calling inappropriate base function of Parameter"); i=false; }
		virtual bool Constrained()				{return false; }
		virtual void SetConstrained(bool i)		{throw ParamExcep("Calling inappropriate base function of Parameter"); i=false; }
		int GetSetting()	{return setting; }
		virtual double GetLowerOfMaxOrUbound()	{throw ParamExcep("Calling inappropriate base function of Parameter"); }
		virtual double GetHigherOfMinOrLbound() {throw ParamExcep("Calling inappropriate base function of Parameter"); }
		virtual std::string GetName()	{std::string c; return c; }

};

class BoundedParameter: public Parameter	{	
	protected :
	double maxv,minv,upbound,lowbound;
	
	public:
	BoundedParameter();
	BoundedParameter(double v,double mn,double mx, int st,double def);
	
	virtual bool HasMin() const;
	virtual void SetMin(double m);
	bool HasMax() const;	
	void SetMax(double m);
	bool HasUpperBound() const;				
	void SetUpperBound(double m)	;
	void SetUpperBound(bool i)		;
	bool HasLowerBound() const;
	void SetLowerBound(double m)	;
	void SetLowerBound(bool i)		;
	bool Constrained()				;
	void SetConstrained(bool i)		;
	void SetCurrent(double v);
	void SetDefault(double m);
	double GetLowerOfMaxOrUbound()	;	 
	double GetHigherOfMinOrLbound() ;
};

class PositiveParameter:public BoundedParameter {
	public :
	PositiveParameter();
	PositiveParameter(double v);
	PositiveParameter(double v,int st,double def);
	PositiveParameter(double v,int st);
	virtual bool HasMin() const {
		return true;
	}
	virtual void SetMin(double m)	{throw ParamExcep("Calling inappropriate base function of PositiveParameter"); m=1; }
	
};	

class NonNegativeParameter:public BoundedParameter	{
	public :
	NonNegativeParameter();
	NonNegativeParameter(double v);
	NonNegativeParameter(double v,int st,double def);
	NonNegativeParameter(double v,int st);
	virtual bool HasMin() const {
		return true;
	}
	virtual void SetMin(double m)	{throw ParamExcep("Calling inappropriate base function of NonNegativeParameter"); m=1; }
	
};	

class FullParameter: public BoundedParameter	{
	std::string name;	
	public :
	FullParameter(double v,double mn,double mx, int st,double def,std::string n);
	std::string GetName() {return name; }
};

class FreqParamGroup	{
	public :
		FreqParamGroup(int n,double *l,int settings);
		FreqParamGroup(int n,Parameter **p);
		~FreqParamGroup();

		const double * const *GetStateFreqsConst() const;

		Parameter *GetParameter(int i);
		void Initialize();
		int GetNParams()	{return nparams; }
		bool IsAMember(Parameter *inp);
		int Fixed();
		int GetSetting();
		
		double GetReparameterized(int n);
		void SetReparameterized(int n,double v);
		double GetReparameterizedMax(int n);
		double GetReparameterizedMin(int n);
		void SetToMinimum(Parameter *d);
		void ForceToSumToOne(double minnonzero);
	private:
		double **GetStateFreqs();

		double sum;
		Parameter **param;
		int nparams;
		double **statefreqs;
		bool owns;
};

inline double **FreqParamGroup::GetStateFreqs()	{
	return statefreqs;
}
inline const double * const * FreqParamGroup::GetStateFreqsConst() const{
	return statefreqs;
}
	
} // namespace bull
#endif
