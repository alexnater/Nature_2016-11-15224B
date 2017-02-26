#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "missing.h"
#include "misc.h"



using namespace std;


MissingData::MissingData(int nbrloci, int nind){
	nloci=nbrloci;
	totsize=nind;
	missdata=new MissdataVector[nloci];
	filearray=new std::string[nloci];
	missingdata=false;
	mmatrix=false;
	}

MissingData::~MissingData(){
	delete [] missdata;
	delete [] filearray;
	if(mmatrix) this->deleteMissMatrix();
	}


int MissingData::ReadMisslist(const char *filelist){
	std::string filename;
	std::ifstream liststream;
	liststream.open( filelist );
	if ( !liststream.is_open() ){
		cout<<"List of missing files "<<filelist<<" could not be opened!\n";
		return 0;
		}
	for (int curlocus=0; curlocus<nloci; ++curlocus){
		if ( liststream.eof() ){ liststream.close(); cout<<"Early end of file in "<<filelist<<"!\n"; return 0; }
		getline(liststream, filename);
		while (filename.size()==0){
			if ( liststream.eof() ){ liststream.close(); cout<<"Early end of file in "<<filelist<<"!\n"; return 0; }
			getline(liststream, filename);
			}
		filearray[curlocus]=filename;
		}
	liststream.close();
	return 1;
	}

int MissingData::ReadMissing(){
	for (int file=0; file<nloci; ++file){
		std::string line;
		ifstream missingfile;
		int ind=0;
		double rangestart, rangeend;
		missdata[file]=MissdataVector(totsize);

		missingfile.open( filearray[file].c_str() );
		if ( !missingfile.is_open() ){
			cout<<"Missing report file "<<filearray[file]<<" could not be opened!\n";
			return 0;
			}

		while ( !missingfile.eof() ){
			getline(missingfile, line);
			if (line.size() != 0){
				istringstream ss(line);
				if ( ss.good() ){
					ss>>ind>>rangestart>>rangeend;
					if (ind>=totsize){
						cout<<"Wrong individual ID "<<ind<<" in missing file: "<<missingfile<<"!\n";
						return 0;
						}
					else {
						missdata[file][ind].push_back(rangestart);
						missdata[file][ind].push_back(rangeend);
						}
					}
				}
			}
//		cout<<"Processed missing report file "<<file+1<<" of "<<nloci<<"."<<endl;
		missingfile.close();
		}
	missingdata=true;
	return 1;
	}

int MissingData::NewInd(int locus, int ind){
	if (locus>=nloci){ return 0; }
	it=missdata[locus][ind].begin();
	end=missdata[locus][ind].end();
	return 1;
	}

bool MissingData::IsMissing(double segpos){
	for (; it != end; it+=2){
		if (segpos<=*it) return false;
		else if (segpos>*it && segpos<=*(it+1) ) return true;
		}
	return false;
	}

bool MissingData::IsMissing(int locus, int ind, double segpos){
	vector<double>::iterator ptr=missdata[locus][ind].begin();
	vector<double>::iterator endptr=missdata[locus][ind].end();
	for (; ptr != endptr; ptr+=2){
		if (segpos<=*ptr) return false;
		else if (segpos>*ptr && segpos<=*(ptr+1) ) return true;
		}
	return false;
	}

bool MissingData::IsMissing(int locus, int ind, int site){
	if (!mmatrix){ cout<<"No missing data matrix available!"<<endl; return 0; }
	else return missmatrix[locus*totsize+ind][site-1];
	}


int MissingData::getValidSites(int nbrsites, int npops, int *popsize, int *minind){
	nsites=nbrsites;
	validsites=ValidSitesVector (nloci);
	for (int locus=0; locus<nloci; ++locus){
		validsites[locus].assign(npops+1, nsites);	// default value if no missing data is provided.
		if (!missingdata) continue;
		for (int site=1; site<=nsites; ++site){
			int allind=0;
			bool accepted=true;
			for (int pop=0; pop<npops; ++pop){
				int covered=0;
				for (int popind=0; popind<popsize[pop]; ++popind, ++allind) if ( !IsMissing(locus, allind, (double)site/(double)nsites) ) ++covered;	// double check correct segpos!!!
				if (covered<minind[pop]){ --validsites[locus][pop]; accepted=false; }	// valid sites per population.
				}
			if (!accepted) --validsites[locus][npops];	// valid sites over all populations.
			}
//		cout<<"Valid sites for locus: "<<locus<<endl;
//		for (int pop=0; pop<=npops; ++pop) cout<<"Population "<<pop<<": "<<validsites[locus][pop]<<"\t";
//		cout<<endl;
		}
	return 1;
	}

int MissingData::createMissMatrix(int nbrsites){
	nsites=nbrsites;
	missmatrix=allocation<bool>(nloci*totsize, nsites);
	for (int locus=0; locus<nloci; ++locus){
		int index=locus*totsize;
		for (int allind=0; allind<totsize; ++allind, ++index){
			for (int site=1; site<=nsites; ++site) if ( !IsMissing(locus, allind, (double)site/(double)nsites) ) missmatrix[index][site-1]=true;
			}
		}
	mmatrix=true;
	return 1;
	}

int MissingData::deleteMissMatrix(){
	if ( !freemem<bool>(missmatrix, nloci*totsize) ){ cout<<"Could not delete matrix of missing data!"<<endl; return 0; }
//	for (int index=0; index<nloci*totsize; ++index) delete [] missmatrix[locus];
//	delete [] missmatrix;
	return 1;
	}




