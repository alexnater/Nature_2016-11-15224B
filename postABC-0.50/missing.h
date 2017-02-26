#ifndef MISSING_H
#define MISSING_H

#include <fstream>


typedef std::vector<std::vector<double> > MissdataVector;
typedef std::vector<std::vector<unsigned int> > ValidSitesVector;

class MissingData {
		char filelist[1000];
		std::string *filearray;
		int totsize;
		int nloci;
		MissdataVector *missdata;
		std::vector<double>::iterator it;
		std::vector<double>::iterator end;

	public:
		MissingData(int nbrloci, int nind);
		~MissingData();
		bool **missmatrix;
		bool missingdata;
		bool mmatrix;
		int nsites;
		ValidSitesVector validsites;
		int ReadMisslist(const char *filelist);
		int ReadMissing();
		int NewInd(int locus, int ind);
		bool IsMissing(double segpos);
		bool IsMissing(int locus, int ind, double segpos);
		bool IsMissing(int locus, int ind, int site);
		int getValidSites(int nbrsites, int npops, int *popsize, int *minind);
		int createMissMatrix(int nbrsites);
		int deleteMissMatrix();

	};

#endif
