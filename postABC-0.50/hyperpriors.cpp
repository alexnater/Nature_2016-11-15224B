#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string.h>
#include "hyperpriors.h"
#include "misc.h"
#include "randlib.h"


using namespace std;


HyperPriors::HyperPriors(int nbrpars, int nbrloci) {
	npars=nbrpars;
	nloci=nbrloci;
	dist=allocation<double>(npars, nloci);
	}

HyperPriors::~HyperPriors() {
	freemem<double>(dist, npars);
	}


void HyperPriors::normtrun(int n, double *array, double mean, double std, double lower, double upper){
	int found=0;
	double cand;
	while (found < n){
		cand=gennor(mean, std);
		if (cand >= lower && cand <= upper) array[found++]=cand;
		}
	}


int HyperPriors::hyperpriors(int Nref, int length, char **params) {

	for (int i=0; i<npars; ++i){
		string type(params[5*i]);
		if (type == "constant"){ double value=atof(params[5*i+1]); for (int j=0; j<nloci; ++j) dist[i][j]=value; }
		else if (type == "beta"){ double scale=atof(params[5*i+1]), a=atof(params[5*i+2]), b=atof(params[5*i+3]); for (int j=0; j<nloci; ++j) dist[i][j]=scale*genbet(a, b); }
		else if (type == "gamma"){ double mean=atof(params[5*i+1]), shape=atof(params[5*i+2]), beta=shape/mean; for (int j=0; j<nloci; ++j) dist[i][j]=gengam(beta, shape); }
		else if (type == "normal"){ double mean=atof(params[5*i+1]), std=atof(params[5*i+2]); for (int j=0; j<nloci; ++j) dist[i][j]=gennor(mean, std); }
		else if (type == "normtrun"){ double mean=atof(params[5*i+1]), std=atof(params[5*i+2]), lower=atof(params[5*i+3]), upper=atof(params[5*i+4]);
			normtrun(nloci, dist[i], mean, std, lower, upper); }
		else { cout<<"Distribution "<<type<<" not supported!"<<endl; exit(0); }
		}

	return 1;
	}


int HyperPriors::printpriors(char *outfile) {
	ofstream file;
	bool usefile=false;
	if (string(outfile) != "stdin" && string(outfile) != "-"){
		usefile=true;
		file.open(outfile);
		if ( !file.is_open() ){
			cout<<"File "<<outfile<<" could not be opened!\n";
			return 0;
			}
		}
	ostream &output=usefile ? file : cout;

	for (int i=0; i<nloci; ++i){
		for (int j=0; j<npars; j++) output<<dist[j][i]<<"\t";
		output<<endl;
		}

	if (usefile) file.close();
	return 1;
	}



