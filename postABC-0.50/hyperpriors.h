#ifndef HYPERPRIORS_H
#define HYPERPRIORS_H


class HyperPriors {

	public:
		int npars, nloci;
		double **dist;

		HyperPriors(int nbrpars, int nbrloci);
		~HyperPriors();
		void normtrun(int n, double *array, double mean, double std, double lower, double upper);
		int hyperpriors(int Nref, int length, char **params);
		int printpriors(char *outfile);
	};


#endif
