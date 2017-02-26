#ifndef INPUT_H
#define INPUT_H

#include <vector>
#include "missing.h"
#include "sumstats.h"


class RawData {

	public:
		double **positions;
		int *nsegs, *popsize;
		char ***haplotypes;
		char ***subhaplotypes;
		int nloci, npops, totsize, subsamplingmode;
		bool usems, subsampling;

		RawData(int maxloci, int npops, int *popsize);
		~RawData();
		int readms_c(char *infile);
		int processms(MissingData &missdata, int nbrloci, bool subsamp, bool ssbefore, int submode, bool haplotize, int *subind, Sumstats &sumstats);
		int subMask(int segsites, int *subind, bool haplotize, short unsigned int **sampledind);
		int subMask(int locus, int segsites, int *subind, bool haplotize, short unsigned int **sampledind);
		int subsample(int popsize, int samplesize, std::vector<int> &selind);

	};

#endif
