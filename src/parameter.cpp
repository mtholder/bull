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
 
#include <cassert>
#include "parameter.hpp"

using namespace bull;

FullParameter::FullParameter(double v,double mn,double mx, int st,double def,std::string n)
	: BoundedParameter(v,mn,mx,st,def)
{name=n;
if (setting&1 && defaultv < minv)	throw ParamExcep(name,"'s default set below its minimum");
if (setting&2 && defaultv>maxv) throw ParamExcep(name,"'s default set above its maximum");
}
PositiveParameter::PositiveParameter(double v,int st,double def)
	: BoundedParameter(v,bull::BULL_SMALL_DBL,0.0,st|par(MIN),def)
{	
	assert(!(setting&par(MAX)));
	assert(v >= 0.0 && def >= 0.0);
}
BoundedParameter::BoundedParameter(double v,double mn,double mx, int st,double def)
	: Parameter(v, st,def)
{	maxv=upbound=mx;
	minv=lowbound=mn;
}
int Parameter::StartWithDefault()	
{return setting&par(DEF);
}
void Parameter::SetStartWithCurrent(bool i) 
{
if (i)	{
		setting&=(~par(ALLSTART));
		setting|=par(CUR);
		}
else	setting&=(~par(CUR));
}
void Parameter::SetDefault(double m)	
{	defaultv=m; 
}

void BoundedParameter::SetCurrent(double m) 
{	if ((setting&par(MIN) && m < minv) || (setting&par(MAX) && m>maxv))
		if (setting&par(MIN) && m < minv)
			if (m<(minv-bull::BULL_SMALL_DBL))
				throw ParamOutOfRangeExcep();
			else
				m=minv;
		else
			if (m>maxv+bull::BULL_SMALL_DBL)
				throw ParamOutOfRangeExcep();
			else
				m=maxv;

	val=m; 
}

void BoundedParameter::SetDefault(double m) 
{ 
defaultv=m; 
if (setting&par(MIN) && defaultv < minv)	throw ParamExcep("Default set below its minimum");
if (setting&par(MAX) && defaultv>maxv)	throw ParamExcep("Default set above its maximum");
}


bool BoundedParameter::HasMin() const {
	return setting&par(MIN);
}

void Parameter::SetStartWithRandom(bool i)	
{
if (i)	{setting&=(~par(ALLSTART));
		setting|=par(RAN);
		}
else	setting&=(~par(RAN));
}

void Parameter::SetStartWithDefault(bool i) 
{	if (i)	{setting&=(~par(ALLSTART));
			setting|=par(DEF);
			}
	else	setting&=(~par(DEF));
}

void Parameter::SetToDefault()	
{	val=defaultv;
}


void Parameter::SetStartWithApproximation(bool i)	
{
if (i)	{setting&=(~par(ALLSTART));
		setting|=par(APPRO);
		}
else	setting&=(~par(APPRO));
}


void Parameter::SetFixed(double m)	
{ 
val=m;
setting= setting|=par(FIX);
}
	

bool BoundedParameter::HasUpperBound()	const {
	return setting&par(UBOUN);
}
void BoundedParameter::SetUpperBound(double m)	
{ 
	upbound=m;
	if (val>upbound)	val=upbound;
	setting= setting|=par(UBOUN);
	if (setting&par(MIN) && upbound < minv) throw ParamExcep("Upper bound set below its minimum");
	if (setting&par(MAX) && upbound>maxv)	throw ParamExcep("Upper bound above its maximum");
}
void BoundedParameter::SetUpperBound(bool i)	
{ 
	assert(!i);
	setting&=(~par(UBOUN));
}

bool BoundedParameter::HasLowerBound() const
{return setting&par(LBOUN);
}
void BoundedParameter::SetLowerBound(double m)	
{ 
	lowbound=m;
	if (val < lowbound) val=lowbound;
	setting= setting|=par(LBOUN);
	if (setting&par(MIN) && lowbound < minv)	throw ParamExcep("Lower bound set below its minimum");
	if (setting&par(MAX) && lowbound>maxv)	throw ParamExcep("Lower bound above its maximum");
}

void BoundedParameter::SetLowerBound(bool i)	
{ 
	assert(!i);
	setting&=(~par(LBOUN));
}


void BoundedParameter::SetConstrained(bool i)	
{
	if (i)	{setting|=par(CONSTR);
			}
	else	setting&=(~par(CONSTR));
}


double BoundedParameter::GetLowerOfMaxOrUbound()	
{	if (par(UBOUN)&setting) return upbound;
	assert(par(MAX));
	return maxv;
}	 

double BoundedParameter::GetHigherOfMinOrLbound()	
{	if (par(LBOUN)&setting) return lowbound;
	assert(par(MIN));
	return minv;

}


Parameter::Parameter(double v,int st,double def)
{	val=v;
	setting=st;
	defaultv=def;
}
void BoundedParameter::SetMin(double m) 
{ 
minv=m; 
if (minv>lowbound)	lowbound=minv;
setting |=par(MIN);
}

bool BoundedParameter::HasMax() const {
	return setting&par(MAX);
}

void BoundedParameter::SetMax(double m) 
{ 
maxv=m; 
if (maxv < upbound) upbound=maxv;
setting |=par(MAX);
}

bool BoundedParameter::Constrained()	
{ return setting&par(CONSTR); 
}
void FreqParamGroup::ForceToSumToOne(double minnonzero)
{	Parameter **temp;
	temp=param;
	sum=0.0;
	for (int i=0; i < nparams; i++)
		{if(param[i]->val>minnonzero)
			sum+=param[i]->val;
		else
			param[i]->val=0.0;
		}
	assert(sum>.999 && sum < 1.001);
	for (int i=0; i < nparams; i++)
		if (param[i]->val>minnonzero)
			param[i]->val/=sum;
	
}	
FreqParamGroup::FreqParamGroup(int n,Parameter **p) 
{
	owns=false;
	statefreqs=new double *[n];
	param=new Parameter *[n];
	nparams=n;
	sum=1.0;
	for (int i=0; i < n; i++)
		{statefreqs[i]=&((*p)->val);
		param[i]=*p++;
		}
}
	

int Parameter::StartWithCurrent()	
{return setting&par(CUR);
}
int Parameter::StartWithRandom()	
{return setting&par(RAN);
}
int Parameter::StartWithApproximation() 
{return setting&par(APPRO) ;
}



FreqParamGroup::~FreqParamGroup()	
{
	if (owns)	for (int i=0; i < nparams; i++)	delete param[i];
	delete []param;
	delete []statefreqs;
}
	

void FreqParamGroup::Initialize()	
{	if ((*param)->StartWithCurrent()) ;
	else if ((*param)->StartWithRandom())
				throw ParamExcep("Random Function to initialize parameters isn't available yet");
			//DirichletDistribution
	else if ((*param)->StartWithApproximation())
				throw ParamExcep("Initial approximation of parameters isn't available yet");
	else if ((*param)->StartWithDefault())
				for (int i=0; i < nparams; i++)
					param[i]->SetToDefault();
	else	throw ParamExcep("No starting value of a parameter has been defined");
}



