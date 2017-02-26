#ifndef TRandomGenerator_H_
#define TRandomGenerator_H_

#include <cmath>
#include <time.h>
#include "TException.h"

class TRandomGenerator{
public:
	long _Idum;
	long usedSeed;

	TRandomGenerator(long addToSeed){
		init(addToSeed);
	};
	TRandomGenerator(long addToSeed, bool seedIsFixed){
		if(!seedIsFixed) init(addToSeed);
		else {
			if(addToSeed<0) addToSeed=-addToSeed;
	        usedSeed=addToSeed;
	        _Idum=-addToSeed;
		}
	};
	TRandomGenerator(){
			init(0);
	};
	~TRandomGenerator(){};
	void init(long addToSeed){
		_Idum=get_randomSeedFromCurrentTime(&addToSeed);
        if(_Idum==0) _Idum=get_randomSeedFromCurrentTime(&_Idum);
        if(_Idum < 0) _Idum=-_Idum;
        if(_Idum > 161803398) _Idum = _Idum % 161803397;
        usedSeed=_Idum;
        _Idum=-_Idum;
	};
	double getRand(){ return ran3(&_Idum);};
	double getRand(double min, double max);
	int getRand(int min, int maxPlusOne);
	long getRand(long min, long maxPlusOne);
	double getNormalRandom (double dMean, double dStdDev);
	long get_randomSeedFromCurrentTime(long* addToSeed);
	float getBiomialRand(float pp, int n);
	float gammln(float xx);
	double getGammaRand(double a, double b);
	double getGammaRand(double a);
	double getGammaRand(int ia);
	double getBetaRandom (double alpha, double beta, double a, double b);
	double getBetaRandom (double alpha, double beta);

private:
	double ran3(long *idum);
};

#endif /* TRandomGenerator_H_ */
