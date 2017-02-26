#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include "input.h"
#include "missing.h"
#include "sumstats.h"
#include "misc.h"
#include "math.h"



using namespace std;

RawData::RawData(int maxloci, int nbrpops, int *ptrpsize){
	usems=true, subsampling=false;
	positions=new double * [maxloci];
	nsegs=new int [maxloci];
	haplotypes=new char ** [maxloci]();
	npops=nbrpops;
	popsize=ptrpsize;
	totsize=sum<int>(npops, popsize);
	nloci=0;
	}

RawData::~RawData(){
	if (usems){
		freemem_c<double>(positions, nloci);
		freemem_c<char>(haplotypes, nloci, totsize);
		}
	else {
		freemem<double>(positions, nloci);
		freemem<char>(haplotypes, nloci, totsize);
		}
	delete [] nsegs;
	if (subsampling && subsamplingmode==0){
		for (int i=0; i < nloci; ++i) delete [] subhaplotypes[i];
		delete [] subhaplotypes;
		}
	}


int RawData::readms_c(char *infile){
	int pop=0, ind=0, allind=0, segsites=0;
	int locus=-1;

	FILE *input;
	if (string(infile) != "stdin" && string(infile) != "-"){
		input=fopen(infile, "r");
		if (input==NULL){
			cout<<"File "<<infile<<" could not be opened!\n";
			return 0;
			}
		}
	else input=stdin;

	char temp[256];
	while(1){
		if ( fscanf(input, " %255s", temp)==EOF ) break;
		else if (strcmp(temp, "segsites:")==0){
			if ( fscanf(input, " %d", &segsites)==EOF ) return 0;
			++locus; ind=0; allind=0; pop=0;
			nsegs[locus]=segsites;
			if (segsites > 0){
				positions[locus]=new double [segsites]();
				haplotypes[locus]=allocation<char>(totsize, segsites+1);
				while(1){
					if ( fscanf(input, " %255s", temp)==EOF ) return 0;
					else if (strcmp(temp, "positions:")==0) break;
					}
				for (int i=0; i<segsites; ++i) if ( fscanf(input, " %lf", positions[locus]+i)==EOF ) return 0;
				char *line;
				while(1){
					line=haplotypes[locus][allind];
					if ( fscanf(input, " %s", line)==EOF ) return 0;
					else {
						if ( strlen(line) != segsites ){
							cout<<"Locus "<<locus<<": length of haplotypes does not match segsites information: "<<strlen(line)<<" versus "<<segsites<<"!"<<endl;
							return 0;
							}
						if ( ind < popsize[pop]-1 ) {++ind; ++allind;}
						else if ( pop < npops-1 ) {ind=0; ++pop; ++allind;}
						else break;
						}
					}
				}
			}
		}
	if (string(infile) != "stdin" && string(infile) != "-") fclose(input);
	usems=false;
	nloci=locus+1;
	return nloci;
	}


int RawData::processms(MissingData &missdata, int nbrloci, bool subsamp, bool ssbefore, int submode, bool haplotize, int *subind, Sumstats &sumstats){
	int subtotsize=sum<int>(npops, subind);
	nloci=nbrloci;
	int length=missdata.nsites;
	subsamplingmode=submode;
	if (subsamp && submode==0){subsampling=true; subhaplotypes=new char **[nloci]; }
	else if (!subsamp) subhaplotypes=haplotypes;

	for (int locus=0; locus<nloci; ++locus){
		short unsigned int **sampledind;
		int segsites=nsegs[locus];
		sumstats.allocate(locus, npops*segsites);
		if (segsites==0) continue;
		int quartile=1;
		for (int seg=0; seg<segsites; ++seg){
			if (positions[locus][seg]>(double)quartile*0.25){ sumstats.quartiles[3*locus+quartile-1]=seg; ++quartile; }
			if (quartile==4) break;
			}
		vector<short unsigned int> indcount(segsites, 0);
		if (subsamp && submode==0){
			subhaplotypes[locus]=new char *[subtotsize];
			sampledind=allocation<short unsigned int>(segsites, totsize);
			subMask(locus, segsites, subind, haplotize, sampledind);
			}
		else if (subsamp && submode==1){
			sampledind=allocation<short unsigned int>(segsites, totsize);
			subMask(segsites, subind, haplotize, sampledind);
			}
		int ind=0, pop=0;
		for (int allind=0; allind<totsize; ++allind){
			char *line=haplotypes[locus][allind];
			if (strlen(line) != segsites ){
				cout<<"Length of haplotype "<<allind<<" at locus "<<locus<<" does not match segsites information: "<<strlen(line)<<" versus "<<segsites<<" !"<<endl;
				if (subsamp) freemem(sampledind, segsites);
				return 0;
				}
			if (missdata.missingdata && !subsamp && !missdata.mmatrix) missdata.NewInd(locus, allind);
			for (int seg=0; seg<segsites; ++seg){
				if (subsamp && !sampledind[seg][allind]) continue;
				bool miss=false;
				if (missdata.mmatrix){
					double site=floor(positions[locus][seg]*length + 0.5);	// round site to nearest integer value.
					miss=missdata.IsMissing(locus, allind, (int)site);
					}
				else if (missdata.missingdata){
					double segpos=positions[locus][seg];
					if (subsamp) miss=missdata.IsMissing(locus, indcount[seg], segpos);
					else miss=missdata.IsMissing(segpos);
					}
				int index=seg*npops+pop;
				if (!miss && line[seg]=='1') ++sumstats.allelemap[locus][index];
				else if (!miss && (line[seg]=='2' || line[seg]=='3') ){ cout<<"Mismatch detected between missing info and genotype state "<<line[seg]<<" at locus "<<locus<<" individual "<<allind<<" segsite "<<seg<<endl; ++sumstats.missing[locus][index]; line[seg]='2'; }
				else if (miss){ ++sumstats.missing[locus][index]; line[seg]='2'; }
				++indcount[seg];
				}
			if ( ind < popsize[pop]-1 ) ++ind;
			else if ( pop < npops-1 ) {ind=0; ++pop;}
			else {
				if (subsamp) freemem(sampledind, segsites);
				if (allind+1 != totsize){ cout<<"Wrong number of samples!"<<endl; return 0; }
				}
			}
//		cout<<"Locus: "<<locus<<" - Segsites: "<<segsites<<endl;
//		for (int seg=0; seg<segsites; ++seg) for (int pop=0; pop<npops; ++pop) cout<<"Segsite: "<<seg<<" Pop: "<<pop<<
//			" Derived: "<<sumstats.allelemap[locus][seg*npops+pop]<<" Missing: "<<sumstats.missing[locus][seg*npops+pop]<<endl;
		}
	return 1;
	}


int RawData::subMask(int segsites, int *subind, bool haplotize, short unsigned int **sampledind){
	int currentsize=0;
	for (int pop=0; pop<npops; pop++){
		int nind=popsize[pop]/2;
		int samplesize;
		if (haplotize){
			samplesize=subind[pop];
			if (samplesize>popsize[pop]/2){
				cout<<"To many individuals selected for subsampling!"<<endl;
				return 0;
				}
			}
		else {
			samplesize=subind[pop]/2;
			if (samplesize>popsize[pop]/2){
				cout<<"To many individuals selected for subsampling!"<<endl;
				return 0;
				}
			}
		int *indices=new int[nind];
		for (int index=0; index<nind; index++) indices[index]=index;

		for (int seg=0; seg<segsites; seg++){
			int counter=nind;
			int sampled=0;
			while (sampled<samplesize){
				int random=rand() % counter;
				int selindex=indices[random];
				if (haplotize){
					int hap=rand() % 2;
					int ind=2*selindex+hap+currentsize;
					sampledind[seg][ind]=1;
					}
				else {
					int ind=2*selindex+currentsize;
					sampledind[seg][ind]=1;
					sampledind[seg][ind+1]=1;
					}
				indices[random]=indices[counter-1];
				indices[counter-1]=selindex;
				sampled++; counter--;
				}
			}
		delete[] indices;
		currentsize+=popsize[pop];
		}
	return 1;
	}

int RawData::subMask(int locus, int segsites, int *subind, bool haplotize, short unsigned int **sampledind){
	int currentsize=0, selind=0;
	for (int pop=0; pop<npops; pop++){
		int nind=popsize[pop]/2;
		int samplesize;
		if (haplotize){
			samplesize=subind[pop];
			if (samplesize>popsize[pop]/2){
				cout<<"To many individuals selected for subsampling!"<<endl;
				return 0;
				}
			}
		else {
			samplesize=subind[pop]/2;
			if (samplesize>popsize[pop]/2){
				cout<<"To many individuals selected for subsampling!"<<endl;
				return 0;
				}
			}
		int *indices=new int[nind];
		for (int index=0; index<nind; index++) indices[index]=index;
		int counter=nind, sampled=0;
		while (sampled<samplesize){
			int random=rand() % counter;
			int selindex=indices[random];
			if (haplotize){
				int hap=rand() % 2;
				int ind=2*selindex+hap+currentsize;
				for (int seg=0; seg<segsites; seg++) sampledind[seg][ind]=1;
				}
			else {
				int ind=2*selindex+currentsize;
				subhaplotypes[locus][selind++]=haplotypes[locus][ind];
				subhaplotypes[locus][selind++]=haplotypes[locus][ind+1];
				for (int seg=0; seg<segsites; seg++){sampledind[seg][ind]=1; sampledind[seg][ind+1]=1;}
				}
			indices[random]=indices[counter-1];
			indices[counter-1]=selindex;
			++sampled; --counter;
			}
		delete[] indices;
		currentsize+=popsize[pop];
		}
	return 1;
	}


int RawData::subsample(int popsize, int samplesize, vector<int> &selind){
//	vector<int> indices(popsize);
	int *indices=new int[popsize];
	for (int index=0; index<popsize; index++) indices[index]=index;
	int sampled=0;
	while (sampled<samplesize){
		int random=rand() % popsize;
		int selindex=indices[random];
		selind[sampled]=selindex;
		indices[random]=indices[popsize-1];
		indices[popsize-1]=selindex;
		sampled++; popsize--;
		}
	delete[] indices;

/*
	double random;
	int processed=0;
	while (sampled<samplesize){
		assert(processed<=nind, "Random subsampling failed!");
		random=(double)rand()/(double)RAND_MAX;
		if ((nind-processed)*random >= samplesize-sampled){
			processed++;
		}
		else {
			sampledind[2*processed]=1;
			sampledind[2*processed+1]=1;
			processed++; sampled++;
			}
		}
*/
	return 1;
	}


