//---------------------------------------------------------------------------
#include "TPrior.h"
//------------------------------------------------------------------------------
 //TPrior
TPrior::TPrior(TRandomGenerator* RandomGenerator, std::string Name, bool IsInt, bool Output){
	randomGenerator=RandomGenerator;
	name=Name;
   isInt=IsInt;
   output=Output;
   curValue=0;
   oldValue=0;
}
void TPrior::makeCurValueInt(){
	curValue=round(curValue);
   /*
	if(curValue-(int)curValue<0.5) curValue=(int)curValue;
   else                           curValue=1+(int)curValue;
   if(curValue<1 && curValue>0) curValue=1;
   */
}
void TPrior::writeValue(const double& value, ofstream& file){
   //attention about scientific notations (integers can't read them)
   if(isInt) file << (long) value;
   else file << value;
}
void TPrior::writeCurValue(ofstream& file){
   writeValue(curValue, file);
}
void TPrior::writeHyperPriorGamma(const double& arg, ofstream& file){
   writeValue(randomGenerator->getGammaRand(arg, arg/curValue), file);
}
void TPrior::writeHyperPriorBeta(ofstream& file){
	if(curValue<0.001) writeValue(0.0, file);
	else {
	   if(curValue>1) throw TException("Hyperprior '"+name+"': mean of the beta distribution > 1!", _FATAL_ERROR);
	   double a=0.5+199*curValue;
	   writeValue(randomGenerator->getBetaRandom(a, a*(1-curValue)/curValue), file);
	}
}
void TPrior::writeHyperPriorNormal(const double& stdev, ofstream& file){
   writeValue(randomGenerator->getNormalRandom(curValue, stdev), file);
}
void TPrior::writeHyperPriorNormalPositive(const double& stdev, ofstream& file){
	float val=randomGenerator->getNormalRandom(curValue, stdev);
	while(val<=0) val=randomGenerator->getNormalRandom(curValue, stdev);
    writeValue(val, file);
}
void TPrior::writeHyperPriorLognormal(const double& stdev, ofstream& file){
	   if(stdev<=0.) throw TException ("Hyperprior '"+name+"': stdev of the log normal distribution <=0!", _FATAL_ERROR);
	   double mean=log(curValue) - 0.5 * log(1 + ( (stdev*stdev)/(curValue*curValue) ));
	   double sigma=sqrt( log( ((stdev*stdev)/(curValue*curValue)) +1 ) );
	   writeValue(exp(randomGenerator->getNormalRandom(mean, sigma)), file);
}
void TPrior::writeHyperPriorLognormalParametersInLog10Scale(const double& stdev, ofstream& file){ //base 10!!!!
	if(stdev<=0.) throw TException ("Hyperprior '"+name+"': stdev of the log normal distribution with parameters in log scale <=0!", _FATAL_ERROR);
	writeValue(pow((double) 10.0, randomGenerator->getNormalRandom(curValue, stdev)), file);
}

void TPrior::saveOldValue(){
   oldValue=curValue;
}
void TPrior::resetOldValue(){
   curValue=oldValue;
}
void TPrior::setCurValue(double newValue){
   curValue=newValue;
   if (isInt) makeCurValueInt();
}

//------------------------------------------------------------------------------
//TCombinedPrior
TCombinedPrior::TCombinedPrior(TRandomGenerator* RandomGenerator, std::string Name, std::string Equation, bool IsInt, bool output, TPriorVector* priorVec):TPrior(RandomGenerator, Name, IsInt, output){
	//replace all - by (-1)*
	Equation=stringReplace('-', "+(-1)*", Equation);
	equation=new TEquation_base(Equation, priorVec);
	//equation->what();
}
//------------------------------------------------------------------------------
//TSimplePiror
TSimplePrior::TSimplePrior(TRandomGenerator* RandomGenerator, std::string Name, double min, double max, bool IsInt, bool output):TPrior(RandomGenerator, Name, IsInt, output){
   lowerLimit=min;
   upperLimit=max;
   if(isInt){
	   lowerLimit-=0.5;
	   upperLimit+=0.5;
   }
   if(lowerLimit>upperLimit) throw TException ("Prior '" + name + "' initialized with min > max!", _FATAL_ERROR);
}
void TSimplePrior::changeCurrentValue(){
   changeCurrentValueWithLimits(lowerLimit, upperLimit);
}
double TSimplePrior::getPriorDensity(){
   return getPriorDensityFromValue(curValue);
}
double TSimplePrior::getOldPriorDensity(){
   return getPriorDensityFromValue(oldValue);
}
void TSimplePrior::changeCurrentValueWithLimits(const double& Min, const double& Max){}
double TSimplePrior::getPriorDensityFromValue(const double& value){
   return 0.0;
}
void TSimplePrior::scale(const double& max, const double& min){
	lowerLimit=(lowerLimit-min)/(max-min);
	upperLimit=(upperLimit-min)/(max-min);
}
//------------------------------------------------------------------------------
//uniform prior
TUniformPrior::TUniformPrior(TRandomGenerator* RandomGenerator, std::string Name, double min, double max, bool IsInt, bool output):TSimplePrior(RandomGenerator, Name, min, max, IsInt, output){}
//standard function with no arguments
void TUniformPrior::changeCurrentValue(){
   changeCurrentValueWithLimits(lowerLimit, upperLimit);
}
double TUniformPrior::getPriorDensity(){
   return getPriorDensityFromValue(curValue);
}
void TUniformPrior::changeCurrentValueWithLimits(const double& Min, const double& Max){
   curValue=randomGenerator->getRand(Min, Max);
   if (isInt) makeCurValueInt();
}
double TUniformPrior::getPriorDensityFromValue(const double& value){
   if(value>upperLimit || value<lowerLimit) return 0.0;
   return (1/(upperLimit-lowerLimit));
}
//------------------------------------------------------------------------------
//fixed prior (always the same value)
TFixedPrior::TFixedPrior(TRandomGenerator* RandomGenerator, std::string Name, double val, bool IsInt, bool output):TSimplePrior(RandomGenerator, Name, val, val, IsInt, output){
	curValue=val;
	if (isInt) makeCurValueInt();
}
//standard function with no arguments
void TFixedPrior::changeCurrentValue(){}
double TFixedPrior::getPriorDensity(){
	//always return 1
	return 1.0;
}
void TFixedPrior::changeCurrentValueWithLimits(const double& Min, const double& Max){}
double TFixedPrior::getPriorDensityFromValue(const double& value){
   if(value==curValue) return 1.0;
   return 0.0;
}
//------------------------------------------------------------------------------
//loguniform prior
TLogUniformPrior::TLogUniformPrior(TRandomGenerator* RandomGenerator, std::string Name, double min, double max, bool IsInt, bool output):TSimplePrior(RandomGenerator, Name, min, max, IsInt, output){
   if( lowerLimit < 1e-15) lowerLimit= 1e-15;
}
void TLogUniformPrior::changeCurrentValue(){
   changeCurrentValueWithLimits(lowerLimit, upperLimit);
}
double TLogUniformPrior::getPriorDensity(){
   return getPriorDensityFromValue(curValue);
}
void TLogUniformPrior::changeCurrentValueWithLimits(const double& Min, const double& Max){
   curValue=exp(randomGenerator->getRand(log(Min), log(Max)));
   if (isInt) makeCurValueInt();
}
double TLogUniformPrior::getPriorDensityFromValue(const double& value){
	if(upperLimit==lowerLimit) return 1.0;
	if(value>upperLimit || value<lowerLimit) return 0.0;
   return 1.0/(value*(log(upperLimit) - log(lowerLimit)));
}
//------------------------------------------------------------------------------
//normal prior
TNormalPrior::TNormalPrior(TRandomGenerator* RandomGenerator, std::string Name, double min, double max, double Mean, double Sigma, bool IsInt, bool output):TSimplePrior(RandomGenerator, Name, min, max, IsInt, output){
   mean=Mean;
   sigma=Sigma;
   if(sigma<=0.) throw TException ("Normal prior '" + name + "' initialized with stdev <=0!", _FATAL_ERROR);
   //now take care of the truncation
   inverseTruncatedArea=0;
   //lower truncation
   if(lowerLimit>0) inverseTruncatedArea+=cumulativeDistributionFunction(lowerLimit);
   //upper truncation
   inverseTruncatedArea+=1-cumulativeDistributionFunction(lowerLimit);
   inverseTruncatedArea=1/inverseTruncatedArea;
}
double TNormalPrior::cumulativeDistributionFunction(double x){
  if(x==0.) return 0.;
  return 0.5*complementaryErrorFunction(-0.707106781186547524*(x-mean)/sigma);
}
double TNormalPrior::complementaryErrorFunction(double x){
   //see book "numerical recipes"
   if(x>=0.) return complementaryErrorCheb(x);
   else return 2.0- complementaryErrorCheb(-x);
}
double TNormalPrior::complementaryErrorCheb(double x){
   //see book "numerical recipes"
   double coef[28]={-1.3026537197817094, 6.4196979235649026e-1,
   1.9476473204185836e-2,-9.561514786808631e-3,-9.46595344482036e-4,
   3.66839497852761e-4,4.2523324806907e-5,-2.0278578112534e-5,
   -1.624290004647e-6,1.303655835580e-6,1.5626441722e-8,-8.5238095915e-8,
   6.529054439e-9,5.059343495e-9,-9.91364156e-10,-2.27365122e-10,
   9.6467911e-11, 2.394038e-12,-6.886027e-12,8.94487e-13, 3.13092e-13,
   -1.12708e-13,3.81e-16,7.106e-15,-1.523e-15,-9.4e-17,1.21e-16,-2.8e-17};
   int j;
   double t, ty, tmp, d=0., dd=0.;
   t=2./(2.+x);
   ty=4.*t-2;
   for(j=27;j>0;--j){
	  tmp=d;
	  d=ty*d-dd+coef[j];
	  dd=tmp;
   }
   return t*exp(-x*x+0.5*(coef[0]+ty*d)-dd);
}
void TNormalPrior::changeCurrentValue(){
   changeCurrentValueWithLimits(lowerLimit, upperLimit);
}
double TNormalPrior::getPriorDensity(){
   return getPriorDensityFromValue(curValue);
}
void TNormalPrior::changeCurrentValueWithLimits(const double& Min, const double& Max){
   do{
	  curValue=randomGenerator->getNormalRandom(mean, sigma);
   }while(curValue<Min || curValue>Max);
   if (isInt) makeCurValueInt();
}
double TNormalPrior::getPriorDensityFromValue(const double& value){
	if(upperLimit==lowerLimit) return 1.0;
	if(value>upperLimit || value<lowerLimit) return 0.0;
   //standardize value
   double y = ( value - mean ) / sigma;
   //no calculate density according to standart formulae, where 0.3989...=1/sqrt(2*PI)
   return 0.398942280401433 * exp ( -0.5 * y * y );
}
void TNormalPrior::scale(const double& max, const double& min){
	lowerLimit=(lowerLimit-min)/(max-min);
	upperLimit=(upperLimit-min)/(max-min);
	mean=(mean-min)/(max-min);
	sigma=sigma/(max-min);
	if(sigma<=0.) throw TException ("Normal prior '" + name + "' scaled to a  stdev <=0!", _FATAL_ERROR);
	//now take care of the truncation
	inverseTruncatedArea=0;
	//lower truncation
	if(lowerLimit>0) inverseTruncatedArea+=cumulativeDistributionFunction(lowerLimit);
	//upper truncation
	inverseTruncatedArea+=1-cumulativeDistributionFunction(lowerLimit);
	inverseTruncatedArea=1/inverseTruncatedArea;
}
//------------------------------------------------------------------------------
//lognormal prior
TLogNormalPrior::TLogNormalPrior(TRandomGenerator* RandomGenerator, std::string Name, double min, double max, double Mean, double Sigma, bool IsInt, bool output):TNormalPrior(RandomGenerator, Name, min, max, Mean, Sigma, IsInt, output){
   if(lowerLimit<0.0) throw TException ("Log normal prior '" + name + "' initialized with min < 0!", _FATAL_ERROR);
   //mean and sigma provided is in x-space, but we need it in the log(x) space!
   meanInX=Mean;
   sigmaInX=Sigma;
   if(Sigma<=0.) throw TException ("Log normal prior '" + name + "' initialized with stdev <=0!", _FATAL_ERROR);
   mean=log(Mean) - 0.5 * log(1 + ( (Sigma*Sigma)/(Mean*Mean) ));
   sigma=sqrt( log( ((Sigma*Sigma)/(Mean*Mean)) +1 ) );
   //now take care of the truncation
   inverseTruncatedArea=0;
   //lower truncation
   if(lowerLimit>0) inverseTruncatedArea+=cumulativeDistributionFunction(log(lowerLimit));
   //upper truncation
   inverseTruncatedArea+=1-cumulativeDistributionFunction(log(lowerLimit));
   inverseTruncatedArea=1/inverseTruncatedArea;
}
void TLogNormalPrior::changeCurrentValue(){
   changeCurrentValueWithLimits(lowerLimit, upperLimit);
}
double TLogNormalPrior::getPriorDensity(){
   return getPriorDensityFromValue(curValue);
}
void TLogNormalPrior::changeCurrentValueWithLimits(const double& Min, const double& Max){
   do{
	  curValue=exp(randomGenerator->getNormalRandom(mean, sigma));
   }while(curValue<Min || curValue>Max);
   if (isInt) makeCurValueInt();
}
double TLogNormalPrior::getPriorDensityFromValue(const double& value){
	if(upperLimit==lowerLimit) return 1.0;
	if(value>upperLimit || value<lowerLimit) return 0.0;
	double y = ( log(value) - sigma ) / mean;
	return inverseTruncatedArea*((1/value)*exp ( -0.5 * y * y ));
}
void TLogNormalPrior::scale(const double& max, const double& min){
	lowerLimit=(lowerLimit-min)/(max-min);
	upperLimit=(upperLimit-min)/(max-min);
	meanInX=(meanInX-min)/(max-min);
	sigmaInX=sigmaInX/(max-min);
    if(lowerLimit<0.) throw TException ("Log normal prior '" + name + "' initialized with min < 0!", _FATAL_ERROR);
    //mean and sigma provided is in x-space, but we need it in the log(x) space!
    if(sigmaInX<=0.) throw TException ("Log normal prior '" + name + "' initialized with stdev <=0!", _FATAL_ERROR);
    mean=log(meanInX) - 0.5 * log(1 + ( (sigmaInX*sigmaInX)/(meanInX*meanInX) ));
    sigma=sqrt( log( ((sigmaInX*sigmaInX)/(meanInX*meanInX)) +1 ) );
    //now take care of the truncation
    inverseTruncatedArea=0;
    //lower truncation
    if(lowerLimit>0) inverseTruncatedArea+=cumulativeDistributionFunction(log(lowerLimit));
    //upper truncation
    inverseTruncatedArea+=1-cumulativeDistributionFunction(log(lowerLimit));
    inverseTruncatedArea=1/inverseTruncatedArea;
}
//------------------------------------------------------------------------------
//TODO: add gamma prior
//gamma prior
/*
TGammaPrior::TGammaPrior(TRandomGenerator* RandomGenerator, std::string Name, double min, double max, bool IsInt):TSimplePrior(RandomGenerator, Name, min, max, IsInt){
   throw TException("Gamma Prior not tested!!!!", _FATAL_ERROR);
}
void TGammaPrior::changeCurrentValue(){
   changeCurrentValueWithLimits(lowerLimit, upperLimit);
}
double TGammaPrior::getPriorDensity(){
   return getPriorDensityFromValue(curValue);
}
void TGammaPrior::changeCurrentValueWithLimits(const double& Min, const double& Max){
   do{
	  curValue=gamma_dev(firstParameter, secondParameter);
   }while(curValue<Min || curValue>Max);
   if (isInt) makeCurValueInt();
}
double TGammaPrior::getPriorDensityFromValue(const double& value){
   if(value>upperLimit || value<lowerLimit) return 0.0;
   return(DFGamma(value, secondParameter, firstParameter));
}
*/

//------------------------------------------------------------------------------
//TEquation
//------------------------------------------------------------------------------
TEquation_base::TEquation_base(std::string & s, TPriorVector* priorVec){
	originalEuqation=s;
	init(s, priorVec, &originalEuqation);
}
TEquation_base::TEquation_base(std::string equat, TPriorVector* priorVec, std::string* Original){
	init(equat, priorVec, Original);
}
void TEquation_base::init(std::string & equat, TPriorVector* priorVec, std::string* Original){
	original=Original;
	myString=equat;
	eraseAllWhiteSpaces(myString);
	equat_object=NULL;
	if(equat.empty()) throw TException("Missing equation!", _FATAL_ERROR);
	if(!generateAppropriateObject(myString, &equat_object, priorVec)) throw TException("Problems reading equation '" + (*original) + "'!", _FATAL_ERROR);
	isConstant=equat_object->isConstant;
}

bool TEquation_base::generateAppropriateObject(std::string equat, TEquation_base** equatObj, TPriorVector* priorVec){
	//variables
	std::string first, second;
	char sign;
	double tmp;
	//find first part of equation
	std::string::size_type pos=equat.find_first_of("+-()");
	if(pos!=std::string::npos && equat.substr(pos,1)==")") throw TException("Unbalanced parenthesis in equation'" + (*original) + "'!", _FATAL_ERROR);

	//no such sign found? -> only contains point operations or no operations
	if(pos==std::string::npos){
		if(stringContainsAny(equat, "*/")){
			//no +, - or ( but still an equation
			first=extractBeforeAnyOf(equat, "*/");
			sign=equat.substr(0,1).c_str()[0];
			equat.erase(0,1);
			if(equat.empty() || first.empty()) throw TException("Problems reading equation '" + (*original) + "'!", _FATAL_ERROR);
			(*equatObj)=new TEquation_operator(first, equat, sign, priorVec, original);
			return true;
		} else {
			//must be a parameter or a value
			TPrior* prior=priorVec->getPriorFromName(equat);
			if(prior!=NULL){
				(*equatObj)=new TEquation_prior(equat, prior, original);
				return true;
			}
			else {
				if((!stringContainsOnly(equat, "0123456789.E")) || (!stringContainsNumbers(equat)))
					throw TException("Problems reading equation '" + (*original) + "': the string '"+equat+"' is unknown!", _FATAL_ERROR);
				tmp=atof(equat.c_str());
				if(tmp==HUGE_VAL || tmp==-HUGE_VAL) throw TException("Problems reading equation '" + (*original) + "': the value '"+equat+"' is out of range!", _FATAL_ERROR);
				(*equatObj)=new TEquation_value(equat, tmp, original);
			}
			return true;
		}
	} else {
		//check if pos is a parenthesis
		while(pos<(equat.size()-1) && equat.substr(pos,1)!="+" && equat.substr(pos,1)!="-" && equat.substr(pos,1)!="*" && equat.substr(pos,1)!="/"){
			if(equat.substr(pos,1)=="("){
				//find closing parenthesis
				int numParenthesis=1;
				std::string::size_type oldPos=pos;
				while(numParenthesis>0){
					pos=equat.find_first_of("()", pos+1);
					if(pos==std::string::npos) throw TException("Unbalanced parenthesis in equation'" + (*original) + "'!", _FATAL_ERROR);
					if(equat.substr(pos,1)=="(") ++numParenthesis;
					else --numParenthesis;
				}
				if(oldPos==0 && pos==(equat.size()-1)){
					//string is surrounded by brackets -> remove them!
					equat.erase(0,1);
					equat.erase(pos-1,1);
					return generateAppropriateObject(equat, equatObj, priorVec);
				} else {
					//search + or - operator after parenthesis
					oldPos=equat.find_first_of("+-()", pos+1);
					if(oldPos==std::string::npos || equat.substr(oldPos,1)=="("){
						//does it contain a * or / before or after the equation?
						//before
						oldPos=equat.find_first_of("*/()");
						if(oldPos==std::string::npos || equat.substr(oldPos,1)=="("){
							//after?
							pos=equat.find_first_of("*/()", pos+1);
							if(pos==std::string::npos  || equat.substr(pos,1)=="(" || equat.substr(pos,1)==")"){
								//throw TException("XXX Missing operator in '" + equat + "' in equation'" + (*original) + "'!", _FATAL_ERROR);
							//else {
								//is probably a function. Get name before first bracket!
								pos=equat.find_first_of("(");
								if(pos==std::string::npos || pos==0)
									throw TException("Problems reading equation '" + (*original) + "'!", _FATAL_ERROR);
								//find closing bracket
								oldPos=pos; int numOpen=1;
								while(numOpen>0){
									oldPos=equat.find_first_of(")(", oldPos+1);
									if(oldPos==std::string::npos) throw TException("Unbalanced parenthesis in equation'" + (*original) + "'!", _FATAL_ERROR);
									if(equat.at(oldPos)=='(') ++numOpen;
									else --numOpen;
								}

								//check if there is something AFTER the function
								if(oldPos!=(equat.size()-1))
									throw TException("Problems reading equation '" + (*original) + "', unexpected characters after closing bracket!", _FATAL_ERROR);
								//TODO: add function with two arguments! -> simply search for comma!
								(*equatObj)=new TEquation_function(equat.substr(0, pos), equat.substr(pos+1, oldPos-pos-1), priorVec, original);
								return true;
							}
						} else {
							if(equat.substr(oldPos,1)==")") throw TException("Unbalanced parenthesis in equation'" + (*original) + "'!", _FATAL_ERROR);
							pos=oldPos;
						}
					} else {
						if(equat.substr(oldPos,1)==")") throw TException("Unbalanced parenthesis in equation'" + (*original) + "'!", _FATAL_ERROR);
						pos=oldPos;
					}
				}
			}
		}

		//split and record sign
		first=equat.substr(0, pos);
		sign=equat.substr(pos,1).c_str()[0];
		second=equat.substr(pos+1);

		if(second.empty()) throw TException("Problems reading equation '" + (*original) + "'!", _FATAL_ERROR);
		if(first.empty()){
			if(sign=='+' || sign=='-') first="0";
			else throw TException("Problems reading equation '" + (*original) + "'!", _FATAL_ERROR);
		}
		TEquation_base* tmpObj=new TEquation_operator(first, second, sign, priorVec, original);
		//(*equatObj)=new TEquation_operator(first, second, sign, priorVec, original);
		if(tmpObj->isConstant){
			(*equatObj)=new TEquation_value(equat, tmpObj->getValue(), original);
			delete tmpObj;
		} else {
			(*equatObj)=tmpObj;
		}
		return true;
	}
	return false;
}

TEquation_value::TEquation_value(std::string s, double Val, std::string* Original){
	myString=s;
	val=Val;
	original=Original;
	isConstant=true;
}
TEquation_prior::TEquation_prior(std::string s, TPrior* Prior, std::string* Original){
	myString=s;
	prior=Prior;
	original=Original;
	isConstant=false;
};

TEquation_operator::TEquation_operator(std::string firstString, std::string secondString, char Sign, TPriorVector* priorVec, std::string* Original){
	myString=firstString+Sign+secondString;
	sign=Sign;
	original=Original;
	first=NULL;
	second=NULL;
	if(sign!='+' && sign!='-' && sign!='*' && sign!='/')
		throw TException("Unknwon operator '", sign, "' in equation '" + (*original) +"'!", _FATAL_ERROR);

	if(!generateAppropriateObject(firstString, &first, priorVec)) throw TException("Problems reading equation '" + (*original) + "'!", _FATAL_ERROR);
	if(!generateAppropriateObject(secondString, &second, priorVec)) throw TException("Problems reading equation '" + (*original) + "'!", _FATAL_ERROR);
	isConstant=first->isConstant && second->isConstant;
};

double TEquation_operator::getValue(){
	double tmp;
	switch(sign){
		case '*': tmp=first->getValue() * second->getValue(); break;
		case '/': tmp=second->getValue();
			  	  if(tmp==0.0){
			  		  throw TException("Problem evaluating equation: division by zero in '" + myString +"'!", _FATAL_ERROR);
			  	  }
			  	  tmp=first->getValue() / tmp;
			  	  break;
		case '+': tmp=first->getValue() + second->getValue(); break;
		case '-': tmp=first->getValue() - second->getValue(); break;
	}
	return tmp;
};
TEquation_function::TEquation_function(std::string function, std::string Equat, TPriorVector* priorVec, std::string* Original){
	myString=function+"("+Equat+")";
	original=Original;
	//find appropriate function
	trimString(function);

	//log and sqrt
	if(function=="exp") pt2Function=&TEquation_function::func_exp;
	else if(function=="pow10") pt2Function=&TEquation_function::func_pow10;
	else if(function=="pow2") pt2Function=&TEquation_function::func_pow2;
	else if(function=="log") pt2Function=&TEquation_function::func_log;
	else if(function=="log10") pt2Function=&TEquation_function::func_log10;
	else if(function=="log2") pt2Function=&TEquation_function::func_log2;
	else if(function=="sqrt") pt2Function=&TEquation_function::func_sqrt;

	//rounding
	else if(function=="abs") pt2Function=&TEquation_function::func_abs;
	else if(function=="floor") pt2Function=&TEquation_function::func_floor;
	else if(function=="ceil") pt2Function=&TEquation_function::func_ceil;
	else throw TException("Unknown function '"+function+"' in equation '" + (*original) +"'!", _FATAL_ERROR);

	//generate equation object of function content
	if(!generateAppropriateObject(Equat, &equat_object, priorVec)) throw TException("Problems reading equation '" + (*original) + "'!", _FATAL_ERROR);
	isConstant=equat_object->isConstant;
}

double TEquation_function::func_exp(double val){
	return exp(val);
}
double TEquation_function::func_pow10(double val){
	return pow(10.0, val);
}
double TEquation_function::func_pow2(double val){
	return pow(2.0, val);
}
double TEquation_function::func_log(double val){
	if(val<=0) throw TException("Problems evaluating equation '" + (*original) + "', can not take logarithm of '", val, "'!", _FATAL_ERROR);
	return log(val);
}
double TEquation_function::func_log10(double val){
	if(val<=0) throw TException("Problems evaluating equation '" + (*original) + "', can not take logarithm of '", val, "'!", _FATAL_ERROR);
	return log10(val);
}
double TEquation_function::func_log2(double val){
	if(val<=0) throw TException("Problems evaluating equation '" + (*original) + "', can not take logarithm of '", val, "'!", _FATAL_ERROR);
	return log2(val);
}

double TEquation_function::func_sqrt(double val){
	if(val<0) throw TException("Problems evaluating equation '" + (*original) + "', can not take square root of '", val, "'!", _FATAL_ERROR);
	return log(val);
}
double TEquation_function::func_abs(double val){
	return fabs(val);
}
double TEquation_function::func_ceil(double val){
	return ceil(val);
}
double TEquation_function::func_floor(double val){
	return floor(val);
}
//------------------------------------------------------------------------------
//TRule
//------------------------------------------------------------------------------
TRule::TRule(std::string firstString, char Sign, std::string secondString, TPriorVector* priorVec){
	sign=Sign;
	if(sign!='<' && sign !='>') throw TException("Cannot initialize rule '" + firstString + " " + Sign + " " + secondString + "': unknown sign '"+Sign+"'!", _FATAL_ERROR);
	first=new TEquation_base(firstString, priorVec);
	second=new TEquation_base(secondString, priorVec);
}
bool TRule::passed(){
	if(sign=='>'){
		if(first->getValue() > second->getValue()) return true;
		else return false;
	} else {
		if(first->getValue() < second->getValue()) return true;
		else return false;
	}
}

//------------------------------------------------------------------------------
//TPriorVector
//------------------------------------------------------------------------------
TPriorVector::TPriorVector(std::string fileName, TLog* gotLogFile, TRandomGenerator* RandomGenerator){
	randomGenerator=RandomGenerator;
	logFile=gotLogFile;
   numPrior=0;
   numSimplePrior=0;
   numCombinedPrior=0;
   mapSimplePrior.empty();
   mapCombinedPrior.empty();
   readPriorsAndRules(fileName);
   _Idum=1L;
}
//---------------------------------------------------------------------------
/** read the *est file and fill the array with the priors */
void TPriorVector::readPriorsAndRules(std::string fileName){
	logFile->listFlush("Reading priors and rules from '" + fileName + "' ...");
	// open file
	ifstream is (fileName.c_str()); // opening the file for reading
	if(!is) throw TException("The .est file '" + fileName + "' could not be opened!", _FATAL_ERROR);

	vector<std::string> simplePriorsNames, combinedPriorsNames;
	vector<std::string>::iterator curName, endName;
	// Reading the file
	std::string my_name, myType, my_equation,my_sign;
	int my_isInt;
	bool my_output;
	double my_min, my_max, my_firstParameter, my_secondParameter;
	std::string::size_type pos;
	int curSection=0; // 0: none, 1: Priors, 2:Rules

	// read line by line
	std::string line;
	vector<std::string> vecLine;
	int lineNum=0;
	while(is.good() && !is.eof()){
	  ++lineNum;
	  //------------------------
	  // go to the prior section
	  //------------------------
	  getline(is, line);
	  line=extractBeforeDoubleSlash(line);   // remove the commentated part
	  trimString(line);
	  if(!line.empty()){
		 if(stringContains(line, "[") && stringContains(line, "]")){
			 if(stringContains(line, "[PARAMETERS]")) { curSection=1; continue;}
			 else if(stringContains(line, "[RULES]")){
				 if(curSection==0) throw TException("Section '[PARAMETERS]]' is missing!", _FATAL_ERROR);
				 curSection=2;
				 continue;
			 }
			 else if(stringContains(line, "[COMPLEX PARAMETERS]")){
				 if(curSection==0) throw TException("Section '[PARAMETERS]]' is missing!", _FATAL_ERROR);
				 curSection=3;
				 continue;
			 } else throw TException("Unknown section '"+line+"'!", _FATAL_ERROR);
		}
		//------------------------
		//switch by section
		//------------------------
		switch (curSection){
		   case 1: //------------------------
			   	   //read the priors
			   	   //------------------------
			   	   fillVectorFromStringWhiteSpaceSkipEmpty(line, vecLine);
			   	   if(vecLine.size()!=5 && vecLine.size()!=6 && vecLine.size()!=8) throw TException("Unexpected number of entries in the .est file '"+fileName+"' on line ", lineNum, "!", _FATAL_ERROR);

			   	   //read if it is int
			   	   my_isInt=stringToInt(vecLine[0]);
			   	   if(my_isInt!=0 && my_isInt !=1) throw TException("Unknown isInt value in the .est file '"+fileName+"' on line ", lineNum, "!", _FATAL_ERROR);
			   	   //read name
			   	   my_name=vecLine[1];
			   	   if(my_name.empty()) throw TException("No name given on line ", lineNum, "!", _FATAL_ERROR);
			   	   //read type
			   	   myType=vecLine[2];
			   	   if(myType!="fixed" && myType!="unif" && myType!="logunif" && myType!="norm" && myType!="lognorm")
			   		   TException("Unknown prior type '" + myType +"' on line ", lineNum, "!", _FATAL_ERROR);
			   	   //read min
			   	   my_min=stringToDouble(vecLine[3]);
			   	   //read last to check for output
			   	   if(vecLine[vecLine.size()-1]=="output") my_output=true;
			   	   else {
			   		   if(vecLine[vecLine.size()-1]!="hide") throw TException("output or hide tag is missing in the .est file '"+fileName+"' on line ", lineNum, "!", _FATAL_ERROR);
			   		   my_output=false;
			   	   }

			   	   //first read fixed prior
			   	   if(myType=="fixed"){
			   		   if(vecLine.size()>5) throw TException("Too many entries for a fixed parameter on line ", lineNum, "!", _FATAL_ERROR);
			   		   mapSimplePrior[my_name]= new TFixedPrior(randomGenerator, my_name, my_min, (bool)my_isInt, my_output);
			   	   } else {
			   		//all priors with 6 entries
			   		if(vecLine.size()<6) throw TException("Too few entries on line ", lineNum, "!", _FATAL_ERROR);
			   		   my_max=stringToDouble(vecLine[4]);
			   		   if(my_min>=my_max) throw TException("Min >= Max on line ", lineNum, "!", _FATAL_ERROR);
			   		   if(myType=="unif"){
			   			   if(vecLine.size()>6) throw TException("Too many entries on line ", lineNum, "!", _FATAL_ERROR);
			   			   mapSimplePrior[my_name]= new TUniformPrior(randomGenerator, my_name, my_min, my_max, (bool)my_isInt, my_output);
			   		   } else if(myType=="logunif"){
			   			   if(vecLine.size()>6) throw TException("Too many entries on line ", lineNum, "!", _FATAL_ERROR);
			   			   mapSimplePrior[my_name]= new TLogUniformPrior(randomGenerator, my_name, my_min, my_max, (bool)my_isInt, my_output);
			   		   } else {
			   			   //all priors with 8 entries
			   			   if(vecLine.size()<8) throw TException("Too few entries on line, keep it short you  ", lineNum, "!", _FATAL_ERROR);
			   			   my_firstParameter=stringToDouble(vecLine[5]);
			   			   my_secondParameter=stringToDouble(vecLine[6]);
			   			   if(myType=="norm"){
			   				   mapSimplePrior[my_name]= new TNormalPrior(randomGenerator, my_name, my_min, my_max, my_firstParameter, my_secondParameter, (bool)my_isInt, my_output);
			   			   } else if(myType=="lognorm"){
			   				   mapSimplePrior[my_name]= new TLogNormalPrior(randomGenerator, my_name, my_min, my_max, my_firstParameter, my_secondParameter, (bool)my_isInt, my_output);
			   			   }
			   		   }
			   	   }
			   	   simplePriorsNames.push_back(my_name);
			   	   break;
		   case	2: //------------------------
			   	   // read the rules
			   	   //------------------------
			   	   pos=line.find_first_of("<>");
			   	   if(pos==std::string::npos) throw TException("Rule '" + line + "' on line ", lineNum, " has no sign (either < or >)!", _FATAL_ERROR);
			   	   rules.push_back(new TRule(line.substr(0,pos), line.substr(pos,1).c_str()[0], line.substr(pos+1), this));
			   	   break;
		   case 3: //------------------------
			   	   // read the complex parameters
			   	   //------------------------
			   	   fillVectorFromStringWhiteSpaceSkipEmpty(line, vecLine);
			   	   if(vecLine.size()<4) throw TException("Too few entries on line ", lineNum, "!", _FATAL_ERROR);
			   	   my_isInt=stringToInt(vecLine[0]);
			   	   if(my_isInt!=0 && my_isInt !=1) throw TException("Unknown isInt value on line ", lineNum, "!", _FATAL_ERROR);
			   	   my_name=vecLine[1];
			   	   if(my_name.empty()) throw TException("No name given on line ", lineNum, "!", _FATAL_ERROR);
			   	   trimString(vecLine[2]);
			   	   if(vecLine[2]!="=") throw TException("Unexcpeted sign (instead of '=') on line ", lineNum, "!", _FATAL_ERROR);
			   	   //read whether or not to output
				   if(vecLine[vecLine.size()-1]=="output") my_output=true;
				   else {
					   if(vecLine[vecLine.size()-1]=="hide")  my_output=false;
					   else throw TException("Missing 'output' or 'hide' tag on line ", lineNum, "!", _FATAL_ERROR);
				   }
				   vecLine.erase(vecLine.end()-1);
			   	   //get equation
			   	   concatenateString(vecLine, my_equation, 3);
			   	   // if the values are ok put save the set as combined paramter
			   	   mapCombinedPrior[my_name]= new TCombinedPrior(randomGenerator, my_name, my_equation, (bool)my_isInt, my_output, this);
			   	   combinedPriorsNames.push_back(my_name);
			   	   break;
			}
	  }
	}

	//create iterators
	curSimplePrior=mapSimplePrior.begin();
	curCombinedPrior=mapCombinedPrior.begin();
	numSimplePrior=simplePriorsNames.size();
	numCombinedPrior=combinedPriorsNames.size();

	//create Arrays of simple and combined Prior pointers for fast reference
	simplePriors = new TSimplePrior*[numSimplePrior];
	combinedPriors = new TCombinedPrior*[numCombinedPrior];

    //keep order of est File: combined priors have to be in the same order as in the est file!!!
	endName=simplePriorsNames.end();
	int i=0;
	for(curName=simplePriorsNames.begin(); curName!=endName; ++curName, ++i){
	   simplePriors[i]=mapSimplePrior[*curName];
	}
	i=0;
	for(curName=combinedPriorsNames.begin(); curName!=combinedPriorsNames.end(); ++curName, ++i)
	   combinedPriors[i]=mapCombinedPrior[*curName];
	is.close();
	logFile->write(" done!");
	logFile->conclude(mapSimplePrior.size(), " parameters");
	logFile->conclude(mapCombinedPrior.size(), " combined parameters");
	logFile->conclude(rules.size(), " rules");
} // readPriorsAndRules
//------------------------------------------------------------------------------
void TPriorVector::writeHeader(ofstream& ofs){
   writeHeaderSimplePriors(ofs);
   writeHeaderCombinedPriors(ofs);
}
void TPriorVector::writeHeaderSimplePriors(ofstream& ofs){
	for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=mapSimplePrior.end(); ++curSimplePrior){
		if(curSimplePrior->second->output) ofs << "\t" << curSimplePrior->second->name;
	}
}
void TPriorVector::writeHeaderCombinedPriors(ofstream& ofs){
	for(curCombinedPrior=mapCombinedPrior.begin(); curCombinedPrior!=mapCombinedPrior.end(); ++curCombinedPrior){
		if(curCombinedPrior->second->output) ofs << "\t" << curCombinedPrior->second->name;
	}
}
//------------------------------------------------------------------------------
void TPriorVector::writeParameters(ofstream& ofs){
   writeParametersSimplePriors(ofs);
   writeParametersCombinedPriors(ofs);
}
void TPriorVector::writeParametersSimplePriors(ofstream& ofs){
	for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=mapSimplePrior.end(); ++curSimplePrior){
	   if(curSimplePrior->second->output){
		   ofs << "\t";
		   curSimplePrior->second->writeCurValue(ofs);
	   }
    }
}
void TPriorVector::writeParametersCombinedPriors(ofstream& ofs){
	for(curCombinedPrior=mapCombinedPrior.begin(); curCombinedPrior!=mapCombinedPrior.end(); ++curCombinedPrior){
		if(curCombinedPrior->second->output){
			ofs << "\t";
			curCombinedPrior->second->writeCurValue(ofs);
		}
    }
}
//------------------------------------------------------------------------------
std::string TPriorVector::getSimplePriorName(int & num){
	return simplePriors[num]->name;
}
TPrior* TPriorVector::getPriorFromName(const std::string& name){
   map<std::string,TSimplePrior*>::iterator myCurSimplePrior=mapSimplePrior.find(name);
   if(myCurSimplePrior!=mapSimplePrior.end()){
	   return myCurSimplePrior->second;
   }
   map<std::string,TCombinedPrior*>::iterator myCurCombinedPrior=mapCombinedPrior.find(name);
   if(myCurCombinedPrior!=mapCombinedPrior.end()) return myCurCombinedPrior->second;
   return NULL;
}
TSimplePrior* TPriorVector::getSimplePriorFromName(const std::string& name){
   map<std::string,TSimplePrior*>::iterator myCurSimplePrior=mapSimplePrior.find(name);
   if(myCurSimplePrior!=mapSimplePrior.end()) return myCurSimplePrior->second;
   return NULL;
}
TCombinedPrior* TPriorVector::getCombinedPriorFromName(const std::string& name){
   map<std::string,TCombinedPrior*>::iterator myCurCombinedPrior=mapCombinedPrior.find(name);
   if(myCurCombinedPrior!=mapCombinedPrior.end()) return myCurCombinedPrior->second;
   return NULL;
}
int TPriorVector::getNumberOfSimplePriorFromName(const std::string& name){
   for(int i=0; i<numSimplePrior;++i) if(simplePriors[i]->name==name) return i;
   return -1;
}
//------------------------------------------------------------------------------
bool TPriorVector::isPriorTag(const std::string& name){
	if(isHyperprior(name)) return true;
	return isPrior(name);
}

bool TPriorVector::isHyperprior(const std::string& name){
	char key=name[0];
	switch (key){
	case '%': case '&': case '$': case '!': case '#': case '@':
		int pos = name.find_last_of(key);
		std::string param = name.substr(pos+1);
		curPrior=getPriorFromName(param);
		if(!curPrior) return false;
		else return true;
		break;
	}
	return false;
}
//------------------------------------------------------------------------------
bool TPriorVector::writeCurValueWithHyperprior(const std::string& name, ofstream& file){
   //check for hyperprior
	char key=name[0];
	std::string param=name;
	switch (key){
			case '%':{ // gamma_deviation
				int pos2 = param.find_first_of('%',1);
				std::string arg = param.substr(1, pos2-1);
				double value=getValue(arg);   // check if it is a prior name
				if(value==_nan) value=atof(arg.c_str());
				curPrior=getPriorFromName(param.substr(pos2+1));
				if(!curPrior) return false;
				curPrior->writeHyperPriorGamma(value, file);
				return true;
                }
			case '&': // beta_deviation
				curPrior=getPriorFromName(param.substr(1));
				if(!curPrior) return false;
				curPrior->writeHyperPriorBeta(file);
				return true;
			case '$': {// normal
				int pos2 = param.find_first_of('$',1);
				std::string arg = param.substr(1, pos2-1);
				double value=getValue(arg);   // check if it is a prior name
				if(value==_nan) value=atof(arg.c_str());
				curPrior=getPriorFromName(param.substr(pos2+1));
				if(!curPrior) return false;
				curPrior->writeHyperPriorNormal(value, file);
				return true;}
			case '!': {// normal truncated at 0 --> always positive!
				int pos2 = param.find_first_of('!',1);
				std::string arg = param.substr(1, pos2-1);
				double value=getValue(arg);   // check if it is a prior name
				if(value==_nan) value=atof(arg.c_str());
				curPrior=getPriorFromName(param.substr(pos2+1));
				if(!curPrior) return false;
				curPrior->writeHyperPriorNormalPositive(value, file);
				return true;}
			case '#': {//lognormal
				int pos2 = param.find_first_of('#',1);
				std::string arg = param.substr(1, pos2-1);
				double value=getValue(arg);   // check if it is a prior name
				if(value==_nan) value=atof(arg.c_str());
				curPrior=getPriorFromName(param.substr(pos2+1));
				if(!curPrior) return false;
				curPrior->writeHyperPriorLognormal(value, file);
				return true;}
			case '@': { //lognorm base 10, but values given in logscale!
				int pos2 = param.find_first_of('@',1);
				std::string arg = param.substr(1, pos2-1);
				double value=getValue(arg);   // check if it is a prior name
				if(value==_nan) value=atof(arg.c_str());
				curPrior=getPriorFromName(param.substr(pos2+1));
				if(!curPrior) return false;
				curPrior->writeHyperPriorLognormalParametersInLog10Scale(value, file);
				return true;}
	}
	return false;
}
//------------------------------------------------------------------------------
bool TPriorVector::writeCurValueToFileFromName(const std::string& name, ofstream& file){
   if(writeCurValueWithHyperprior(name, file)) return true;
   curPrior=getPriorFromName(name);
   if(!curPrior) return false;
   curPrior->writeCurValue(file);
   return true;
}
//------------------------------------------------------------------------------
void TPriorVector::getNewValues(){
	bool checkRules=false;
	while(!checkRules){
		// change all simple prior values
		for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=mapSimplePrior.end(); ++curSimplePrior)
			   curSimplePrior->second->changeCurrentValue();
		//check the rules if the rules do not fit change ALL parameters until the rules fit
		checkRules=true;
		for(vector<TRule*>::iterator curRule = rules.begin(); curRule!=rules.end(); ++curRule){
			if(!(*curRule)->passed()) checkRules=false;
		}
	}
   // calc all combined parameters
   updateCombinedParameters();
}
//------------------------------------------------------------------------------
//functions to update parameters via MCMC
void TPriorVector::getNewValuesMcmc(){
	//first get new values
	double* oldValues=new double[numSimplePrior];
	for(int i=0;i<numSimplePrior;++i) oldValues[i]=simplePriors[i]->curValue;
	bool rulesok=false;
	while(!rulesok){
		for(int i=0;i<numSimplePrior;++i){

			double r=randomGenerator->getRand(0.0, simplePriors[i]->mcmcStep);
			double s=simplePriors[i]->mcmcStep/2.0;
			simplePriors[i]->curValue=oldValues[i]+r-s;
			//simplePriors[i]->curValue=oldValues[i]+randomGenerator->getRand(0.0, simplePriors[i]->mcmcStep) - simplePriors[i]->mcmcStep/2.0;
			//reflect, if necessary
			if(simplePriors[i]->curValue < simplePriors[i]->lowerLimit) simplePriors[i]->curValue=simplePriors[i]->lowerLimit+(simplePriors[i]->lowerLimit-simplePriors[i]->curValue);
			if(simplePriors[i]->curValue > simplePriors[i]->upperLimit) simplePriors[i]->curValue=simplePriors[i]->upperLimit+(simplePriors[i]->upperLimit-simplePriors[i]->curValue);
		}
		//check the rules if the rules do not fit change ALL parameters until the rules fit
		rulesok=true;
		for(vector<TRule*>::iterator curRule = rules.begin(); curRule!=rules.end(); ++curRule){
			if(!(*curRule)->passed()) rulesok=false;
		}
	}
	// calc all combined parameters
	updateCombinedParameters();
	delete[] oldValues;
}
void TPriorVector::getNewValuesMcmcUpdateOnePriorOnly(){
   //select prior
   int r=numSimplePrior*randomGenerator->getRand()-0.5;
   getNewValuesMcmc(simplePriors[r]);
   // calc all combined parameters
   updateCombinedParameters();
}
void TPriorVector::getNewValuesMcmc(const int& priorNumber){
	  getNewValuesMcmc(simplePriors[priorNumber]);
	   // calc all combined parameters
	   updateCombinedParameters();
}
void TPriorVector::getNewValuesMcmc(const int& priorNumber, double & mcmcStep){
	  getNewValuesMcmc(simplePriors[priorNumber], mcmcStep);
	   // calc all combined parameters
	   updateCombinedParameters();
}
void TPriorVector::getNewValuesMcmc(TSimplePrior* thisSimplePrior){
	getNewValuesMcmc(thisSimplePrior, thisSimplePrior->mcmcStep);
}

void TPriorVector::getNewValuesMcmc(TSimplePrior* thisSimplePrior, double & mcmcStep){
   bool rulesok=false;
   double oldValue=thisSimplePrior->curValue;
   while(!rulesok){
   		for(int i=0;i<numSimplePrior;++i){
   			thisSimplePrior->curValue=oldValue+randomGenerator->getRand(0.0, mcmcStep) - mcmcStep/2;
   			//reflect, if necessary
   			if(thisSimplePrior->curValue < thisSimplePrior->lowerLimit) thisSimplePrior->curValue=thisSimplePrior->lowerLimit+(thisSimplePrior->lowerLimit-thisSimplePrior->curValue);
   			if(thisSimplePrior->curValue > thisSimplePrior->upperLimit) thisSimplePrior->curValue=thisSimplePrior->upperLimit+(thisSimplePrior->upperLimit-thisSimplePrior->curValue);
   		}
   		//check the rules if the rules do not fit change ALL parameters until the rules fit
   		rulesok=true;
   		for(vector<TRule*>::iterator curRule = rules.begin(); curRule!=rules.end(); ++curRule){
   			if(!(*curRule)->passed()) rulesok=false;
   		}
   }
}
//------------------------------------------------------------------------------
//functions to update parameters via PMC
bool TPriorVector::getNewValuesPMC(double* newParams){
	//save new params into prior objects...
	for(int i=0; i<numSimplePrior; ++i) simplePriors[i]->setCurValue(newParams[i]);

	for(vector<TRule*>::iterator curRule = rules.begin(); curRule!=rules.end(); ++curRule){
		if(!(*curRule)->passed()) return false;
	}

   // calc all combined parameters
   updateCombinedParameters();
   return true;
}

//------------------------------------------------------------------------------
void TPriorVector::resetOldValues(){
	for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=mapSimplePrior.end(); ++curSimplePrior)
	   curSimplePrior->second->resetOldValue();
	for(curCombinedPrior=mapCombinedPrior.begin(); curCombinedPrior!=mapCombinedPrior.end(); ++curCombinedPrior)
	   curCombinedPrior->second->resetOldValue();
};
//------------------------------------------------------------------------------
void TPriorVector::saveOldValues(){
   for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=mapSimplePrior.end(); ++curSimplePrior)
	   curSimplePrior->second->saveOldValue();
   for(curCombinedPrior=mapCombinedPrior.begin(); curCombinedPrior!=mapCombinedPrior.end(); ++curCombinedPrior)
	   curCombinedPrior->second->saveOldValue();
};
//------------------------------------------------------------------------------
void TPriorVector::updateCombinedParameters(){
	//same order as in the est file!!!
	for(int i=0; i<numCombinedPrior; ++i)
		combinedPriors[i]->update();
};
//------------------------------------------------------------------------------
double TPriorVector::getPriorDensity(){
   double density=1.0;
   for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=mapSimplePrior.end(); ++curSimplePrior)
	  density*=curSimplePrior->second->getPriorDensity();
   return density;
};
double TPriorVector::getOldPriorDensity(){
   double density=1.0;
   for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=mapSimplePrior.end(); ++curSimplePrior)
	  density*=curSimplePrior->second->getOldPriorDensity();
   return density;
};
double TPriorVector::getPriorDensity(double* values){
   double density=1.0;
   for(int i=0; i<numSimplePrior; ++i){
	  density*=simplePriors[i]->getPriorDensityFromValue(values[i]);
   }
   return density;
};
//------------------------------------------------------------------------------
double TPriorVector::getValue(const std::string& name){
   curPrior=getPriorFromName(name);
   if(!curPrior) return _nan;
   return curPrior->curValue;
}
//------------------------------------------------------------------------------
bool TPriorVector::isPrior(const std::string& name){
   curPrior=getPriorFromName(name);
   if(!curPrior) return false;
   return true;
}
//------------------------------------------------------------------------------
void TPriorVector::setSimplePriorMCMCStep(const int& priorNumber, const double& prop){
   simplePriors[priorNumber]->mcmcStep=prop*(simplePriors[priorNumber]->upperLimit-simplePriors[priorNumber]->lowerLimit);
}
void TPriorVector::setSimplePriorValue(const int& priorNumber, const double& val){
	simplePriors[priorNumber]->setCurValue(val);
}
void TPriorVector::setSimplePriorValue(const std::string& name, const double& val){
	curPrior=getPriorFromName(name);
	curPrior->setCurValue(val);
}




