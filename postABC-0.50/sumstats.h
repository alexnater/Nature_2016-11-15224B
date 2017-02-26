#ifndef SUMSTATS_H
#define SUMSTATS_H

#include <vector>
#include <string>
#include <map>
#include "missing.h"
#include "TPrior.h"

typedef std::vector<std::vector<short unsigned int> > GroupVector; //groupvector[group][pop]
// typedef std::vector<std::vector<unsigned int> > ValidsitesVector; //validsitesvector[group][locus]
typedef std::vector<std::vector<bool> > ValidSegsVector; //validSegsVector[locus][seg]
typedef std::vector<std::vector<double> > SumstatVector; //sumstats[sumstat][locus*pop]
typedef std::vector<double> OneSumstatVector; //sumstats[locus*pop]
typedef std::map<std::string, OneSumstatVector> SumstatMap; //sumstats[sumstatname]->[locus*pop]
typedef std::vector<std::vector<double> > LocusSummary; //locsummary[sumstat][2*pop]
typedef std::map<std::string, std::vector<double> > LocusSummaryMap; //locsummary[sumstatname][2*pop]


class Sumstats {

	public:
		short unsigned int **allelemap, **missing, *quartiles;
		GroupVector groups;
		ValidSitesVector validsites;
		ValidSegsVector validsegs;
		SumstatVector sumstats;
		SumstatMap sumstatmap;
		LocusSummary locsummary;
		LocusSummaryMap locsummarymap;
		std::vector<std::string> sumstatnames;
		int nloci;

		Sumstats(int nbrloci);
		~Sumstats();
		int readgroupinfo(char *infile);
		int allocate(int locus, int size);
		int getValidSegs(int *nsegs, int npops, int *popsize, int *minind);
		int meansd(std::string name, int period, std::vector<double>::iterator start, std::vector<double>::iterator end);
		int moments(std::string name, int period, std::vector<double>::iterator start, std::vector<double>::iterator end);
		int calculateSegs(int *nsegs, int npops, int *popsize, ValidSitesVector &gotValidsites, int scaling);
		int calcDiversity(int *nsegs, int npops, int *popsize, int *minind, ValidSitesVector &gotValidsites, int scaling);
		int tajD(int npops, int *popsize, bool teststat=false);
		double denominator(int n, int segsites);
		int ZnS(int *nsegs, int npops, int *popsize, int *minind, char ***haplotypes);
		int ZnS_Quartiles(int *nsegs, int npops, int *popsize, int *minind, char ***haplotypes);
		double r2(std::vector<std::pair<int, int> >::iterator start, std::vector<std::pair<int, int> >::iterator end);
		int sfs(bool folded, int *nsegs, int npops, int *popsize, int *minind);
		int printSumstats(char *outfile);
		int printSumstatsPriors(char *outfile, int nsim, TPriorVector *priors, bool printHeader, bool append);
	};
 
#endif
