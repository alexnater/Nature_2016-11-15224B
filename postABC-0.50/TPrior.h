//---------------------------------------------------------------------------
#ifndef TPriorH
#define TPriorH
#include "TException.h"
#include "TRandomGenerator.h"
#include "stringFunctions.h"
#include "TLog.h"
#include <vector>
#include <map>
//TODO replace with a better solution. For instance, pass by reference and return false?
const double _nan=-999999;

//------------------------------------------------------------------------------
//TPrior
//------------------------------------------------------------------------------
class TPrior{
   public:
	  std::string name;
	  bool isInt;
	  bool output;
	  double curValue, oldValue;
	  TRandomGenerator* randomGenerator;

	  TPrior(TRandomGenerator* RandomGenerator, std::string Name, bool isInt, bool Output);
	  ~TPrior(){}
	  void makeCurValueInt();
	  void writeValue(const double& value, ofstream& file);
	  void writeCurValue(ofstream& file);
	  void writeHyperPriorGamma(const double& arg, ofstream& file);
	  void writeHyperPriorBeta(ofstream& file);
	  void writeHyperPriorNormal(const double& stdev, ofstream& file);
	  void writeHyperPriorNormalPositive(const double& stdev, ofstream& file);
	  void writeHyperPriorLognormal(const double& stdev, ofstream& file);
	  void writeHyperPriorLognormalParametersInLog10Scale(const double& stdev, ofstream& file); // base 10!!!
	  void saveOldValue();
	  void resetOldValue();
	  void setCurValue(double newValue);
};

//------------------------------------------------------------------------------
//TSimplePrior
//------------------------------------------------------------------------------
class TSimplePrior:public TPrior{
   public:
	  double mcmcStep;
	  double upperLimit, lowerLimit;

	  TSimplePrior(TRandomGenerator* RandomGenerator, std::string Name, double min, double max, bool IsInt, bool output);
	  virtual ~TSimplePrior(){};
	  virtual void changeCurrentValue();
	  virtual double getPriorDensity();
	  virtual double getOldPriorDensity();
	  virtual void changeCurrentValueWithLimits(const double& Min, const double& Max);
	  virtual double getPriorDensityFromValue(const double& value);
	  virtual void scale(const double& max, const double& min);
};
//------------------------------------------------------------------------------
class TUniformPrior:public TSimplePrior{
   public:
	  TUniformPrior(TRandomGenerator* RandomGenerator, std::string Name, double min, double max, bool IsInt, bool output);
	   void changeCurrentValueWithLimits(const double& Min, const double& Max);
	   void changeCurrentValue();
	   double getPriorDensity();
	   double getPriorDensity(const double& value);
	   double getPriorDensityFromValue(const double& value);
};
class TFixedPrior:public TSimplePrior{
   public:
	TFixedPrior(TRandomGenerator* RandomGenerator, std::string Name, double val, bool IsInt, bool output);
	   void changeCurrentValueWithLimits(const double& Min, const double& Max);
	   void changeCurrentValue();
	   double getPriorDensity();
	   double getPriorDensity(const double& value);
	   double getPriorDensityFromValue(const double& value);
};
class TLogUniformPrior:public TSimplePrior{
   public:
      //double inverseTruncatedArea;
	  TLogUniformPrior(TRandomGenerator* RandomGenerator, std::string Name, double min, double max, bool IsInt, bool output);
	  void changeCurrentValue();
	  double getPriorDensity();
	  void changeCurrentValueWithLimits(const double& Min, const double& Max);
	  double getPriorDensity(const double& value);
	  double getPriorDensityFromValue(const double& value);
};
class TNormalPrior:public TSimplePrior{
   public:
   double mean, sigma, inverseTruncatedArea;
	  TNormalPrior(TRandomGenerator* RandomGenerator, std::string Name, double min, double max, double Mean, double Sigma, bool IsInt, bool output);
	  void changeCurrentValue();
	  double cumulativeDistributionFunction(double x);
	  double complementaryErrorFunction(double x);
	  double complementaryErrorCheb(double x);
	  double getPriorDensity();
	  void changeCurrentValueWithLimits(const double& Min, const double& Max);
	  double getPriorDensity(const double& value);
	  double getPriorDensityFromValue(const double& value);
	  void scale(const double& min, const double& max);
};
class TLogNormalPrior:public TNormalPrior{
   public:
	  double meanInX, sigmaInX;
	  TLogNormalPrior(TRandomGenerator* RandomGenerator, std::string Name, double min, double max, double Mean, double Sigma, bool IsInt, bool output);
	  void changeCurrentValue();
	  double getPriorDensity();
	  void changeCurrentValueWithLimits(const double& Min, const double& Max);
	  double complementaryErrorFunction(double x);
	  double complementaryErrorCheb(double x);
	  double getPriorDensity(const double& value);
	  double getPriorDensityFromValue(const double& value);
	  void scale(const double& min, const double& max);
};
class TGammaPrior:public TSimplePrior{
   public:
	  TGammaPrior(TRandomGenerator* RandomGenerator, std::string Name, double FirstParameter, double SecondParameter, bool IsInt, bool output);
	  void changeCurrentValue();
	  double getPriorDensity();
	  void changeCurrentValueWithLimits(const double& Min, const double& Max);
	  double getPriorDensity(const double& value);
	  double getPriorDensityFromValue(const double& value);
};
//---------------------------------------------------------------------------
//TEquation
//---------------------------------------------------------------------------
class TPriorVector;
class TEquation_base{
private:
	std::string originalEuqation;
public:
	TEquation_base* equat_object;
	std::string* original;
	std::string myString;
	bool isConstant;
	TEquation_base(){equat_object=NULL;isConstant=false;};
	TEquation_base(std::string & s, TPriorVector* priorVec);
	TEquation_base(std::string s, TPriorVector* priorVec, std::string* Original);
	virtual ~TEquation_base(){if(equat_object!=NULL){delete equat_object;}};
	void init(std::string & equat, TPriorVector* priorVec, std::string* Original);
	bool generateAppropriateObject(std::string s, TEquation_base** equatObj, TPriorVector* priorVec);
	virtual double getValue(){return equat_object->getValue();};
	virtual void what(){cout << "TEquation_base from '" << myString << "'"; if(isConstant) cout << "(constant) "; cout << " with:" << endl; equat_object->what();};
	bool replaceOperatorsOfValuesWithValue();

};

class TEquation_value: public TEquation_base{
public:
	double val;
	TEquation_value(std::string s, double Val, std::string* Original);
	~TEquation_value(){};
	double getValue(){return val;};
	virtual void what(){cout << "TEquation_value from '" << myString << "'"; if(isConstant) cout << "(constant)"; cout << endl; };
};

class TEquation_prior: public TEquation_base{
public:
	TPrior* prior;
	TEquation_prior(std::string s, TPrior* Prior, std::string* Original);
	~TEquation_prior(){};
	 double getValue(){return prior->curValue;};
	 virtual void what(){cout << "TEquation_prior from '" << myString << "'"; if(isConstant) cout << "(constant)"; cout << endl; };
};

class TEquation_function: public TEquation_base{
public:
	double (TEquation_function::*pt2Function)(double);

	TEquation_function(std::string function, std::string Equat, TPriorVector* priorVec, std::string* Original);
	~TEquation_function(){};
	//implemented functions (allows for better error handling)
	double func_exp(double val);
	double func_pow10(double val);
	double func_pow2(double val);
	double func_log(double val);
	double func_log10(double val);
	double func_log2(double val);
	double func_sqrt(double val);
	double func_abs(double val);
	double func_ceil(double val);
	double func_floor(double val);

	double getValue(){
		return (*this.*pt2Function) (equat_object->getValue());
	};
	virtual void what(){cout << "TEquation_function from '" << myString << "' "; if(isConstant) cout << "(constant) "; cout << " with:" << endl;
		equat_object->what();
	};
};

class TEquation_operator: public TEquation_base{
public:
	TEquation_base* first;
	TEquation_base* second;
	char sign;

	TEquation_operator(std::string firstString, std::string secondString, char Sign, TPriorVector* priorVec, std::string* Original);
	~TEquation_operator(){
		if(first!=NULL) delete first;
		if(second!=NULL) delete second;
	};
	double getValue();
	virtual void what(){
		cout << "TEquation_operator from '" << myString << "' "; if(isConstant) cout << "(constant) "; cout << "with:" << endl;
		first->what();
		second->what();
	};
};

//------------------------------------------------------------------------------
//TCombinedPrior
//------------------------------------------------------------------------------
class TCombinedPrior: public TPrior{
   public:
		TEquation_base* equation;
		TCombinedPrior(TRandomGenerator* RandomGenerator, std::string Name, std::string Equation, bool IsInt, bool output, TPriorVector* priorVec);
		~TCombinedPrior(){
			delete equation;
		}
		void update(){
			curValue=equation->getValue();
		};
};


//---------------------------------------------------------------------------
//TRule
//---------------------------------------------------------------------------
class TRule{
   public:
	  TEquation_base* first;
	  TEquation_base* second;
	  char sign;

	  TRule(std::string firstString, char Sign, std::string secondString, TPriorVector* priorVec);
	  ~TRule(){
		  if(first) delete first;
		  if(second) delete second;
	  }
	  bool passed();
};
//---------------------------------------------------------------------------
//TPriorVector
//---------------------------------------------------------------------------
//TODO: make sub class using TEstimator to initialize and use it for posterior predictive simulations
class TPriorVector{
   //friend istream& operator >> (istream& is, TPriorVector& S);
	//friend ostream& operator << (ostream& os, const TPriorVector& S);
	public:
	    TRandomGenerator* randomGenerator;
		map<std::string,TSimplePrior*> mapSimplePrior;
		map<std::string,TCombinedPrior*> mapCombinedPrior;
		map<std::string,TSimplePrior*>::iterator curSimplePrior;
		map<std::string,TCombinedPrior*>::iterator curCombinedPrior;
		vector<TRule*> rules;
		int numPrior;
		int numSimplePrior;
		int numCombinedPrior;
		TSimplePrior** simplePriors;
		TCombinedPrior** combinedPriors;
		TPrior* curPrior;
		long _Idum;
		TLog* logFile;

		TPriorVector(std::string fileName, TLog* gotLogFile, TRandomGenerator* randomGenerator);
		~TPriorVector(){
			for(vector<TRule*>::iterator it=rules.begin(); it!=rules.end(); ++it)
				delete (*it);
		   for(curSimplePrior=mapSimplePrior.begin(); curSimplePrior!=mapSimplePrior.end(); ++curSimplePrior){
			   delete curSimplePrior->second;
		   }
		   delete[] simplePriors;
		   for(curCombinedPrior=mapCombinedPrior.begin(); curCombinedPrior!=mapCombinedPrior.end(); ++curCombinedPrior){
			   delete curCombinedPrior->second;
		   }
		   delete[] combinedPriors;
		}
		std::string getSimplePriorName(int & num);
		void getNewValues();
		void resetOldValues();
		void saveOldValues();
		void readPriorsAndRules(std::string fileName);
		void updateCombinedParameters();
		void writeHeader(ofstream& ofs);
		void writeHeaderSimplePriors(ofstream& ofs);
		void writeHeaderCombinedPriors(ofstream& ofs);
		void writeParameters(ofstream& ofs);
		void writeParametersSimplePriors(ofstream& ofs);
		void writeParametersCombinedPriors(ofstream& ofs);
		void getNewValuesMcmcUpdateOnePriorOnly();
		void getNewValuesMcmc(const int& priorNumber);
		void getNewValuesMcmc(const int& priorNumber, double & mcmcStep);
		void getNewValuesMcmc();
		void getNewValuesMcmc(TSimplePrior* thisSimplePrior);
		void getNewValuesMcmc(TSimplePrior* thisSimplePrior, double & mcmcStep);
		bool getNewValuesPMC(double* newParams);
		double getPriorDensity();
		double getOldPriorDensity();
		double getPriorDensity(double* values);
		double getValue(const std::string& Name);
		bool isPrior(const std::string& name);

		double calcEquation(const std::string & equation);
		double calcEquation(std::string equat, const std::string & original);
		bool isPriorTag(const std::string& name);
		bool isHyperprior(const std::string& name);
		bool writeCurValueToFileFromName(const std::string& name, ofstream& file);
	    bool writeCurValueWithHyperprior(const std::string& name, ofstream& file);
		TPrior* getPriorFromName(const std::string& name);
		TSimplePrior* getSimplePriorFromName(const std::string& name);
		TCombinedPrior* getCombinedPriorFromName(const std::string& name);
		int getNumberOfSimplePriorFromName(const std::string& name);
		void setSimplePriorValue(const int& priorNumber, const double& value);
		void setSimplePriorValue(const std::string& name, const double& val);
		void setSimplePriorMCMCStep(const int& priorNumber, const double& prop);
};



#endif
