#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <new>
#include "sumstats.h"
#include "TPrior.h"
#include "misc.h"



using namespace std;


Sumstats::Sumstats(int nbrloci){
	nloci=nbrloci;
	allelemap=new short unsigned int * [nloci];
	missing=new short unsigned int * [nloci];
	quartiles=new short unsigned int [3*nloci]();
	}

Sumstats::~Sumstats(){
	freemem<short unsigned int>(allelemap, nloci);
	freemem<short unsigned int>(missing, nloci);
	}


int Sumstats::readgroupinfo(char *infile){	// under construction!
	std::ifstream groupinfo;
	int groupindex;
	std::string popindices;
	std::string nvalid;

	groupinfo.open( infile );
	if ( !groupinfo.is_open() ){
		cout<<"File with number of valid sites "<<infile<<" could not be opened!\n";
		return 0;
		}
	while ( groupinfo>>groupindex>>popindices>>nvalid ){
		groups.push_back( vector<short unsigned int>() );
		std::istringstream ss1(popindices);
		while ( !ss1.eof() ){
			std::string field;
			getline( ss1, field, ',' );
			groups.back().push_back( atoi( field.c_str() ) );
			}
		validsites.push_back( vector<unsigned int>(nloci, 0) );
		std::istringstream ss2(nvalid);
		while ( !ss2.eof() ){
			std::string field;
			getline( ss2, field, ',' );
			validsites.back().push_back( atoi( field.c_str() ) );
			}
		}
	groupinfo.close();
	return 1;
	}


int Sumstats::allocate(int locus, int size){
	try {
		allelemap[locus]=new short unsigned int [size]();
		missing[locus]=new short unsigned int [size]();
		}
	catch (std::bad_alloc& ba){
		std::cerr << "bad_alloc caught: " << ba.what() << '\n';
		}
	return 1;
	}


int Sumstats::meansd(std::string name, int npops, vector<double>::iterator start, vector<double>::iterator end){
	locsummarymap[name]=std::vector<double>(2*npops, 0);
	std::vector<double> &itpopsum=locsummarymap[name];

	if (start==end){
		cout<<"Empty vector of summary statistics!"<<endl;
		return 0;
		}

	for (int pop=0; pop<npops; ++pop){
		double rt=0, sum2=0, am, sd;
		int n=0;
		for (vector<double>::iterator it=start+pop; it != end+pop; it+=npops){ if (isfinite(*it)){ rt+=*it; ++n; } }
		am=rt/(double)n;

		for (vector<double>::iterator it=start+pop; it != end+pop; it+=npops){ if (isfinite(*it)) sum2+=pow(((*it)-am), (double)2); }
		sd=sqrt(sum2/(double)n);

		itpopsum[2*pop]=am;
		itpopsum[2*pop+1]=sd;
		}

	return 1;
	}


int Sumstats::moments(std::string name, int npops, vector<double>::iterator start, vector<double>::iterator end){
	locsummarymap[name]=std::vector<double>(4*npops, 0);
	std::vector<double> &itpopsum=locsummarymap[name];

	if (start==end){
		cout<<"Empty vector of summary statistics!"<<endl;
		return 0;
		}

	for (int pop=0; pop<npops; ++pop){
		double rt=0, sum2=0, sum3=0, sum4=0;
		int n=0;
		for (vector<double>::iterator it=start+pop; it != end+pop; it+=npops){ if ( isfinite(*it) ){ rt+=*it; ++n; } }
		double am=rt/(double)n;

		for (vector<double>::iterator it=start+pop; it != end+pop; it+=npops){
			if ( isfinite(*it) ){
				double diff=(*it)-am;
				double diff2=diff*diff;
				sum2+=diff2;
				sum3+=diff2*diff;
				sum4+=diff2*diff2;
				}
			}
		double var=sum2/(double)n;
		double sd=sqrt(var);
		double skew=(var>0) ? (sum3/(double)n)/pow(sd, (double)3) : 0;
		double kurt=(var>0) ? ((sum4/(double)n)/pow(var, (double)2))-3 : 0;

		itpopsum[4*pop]=am;
		itpopsum[4*pop+1]=var;
		itpopsum[4*pop+2]=skew;
		itpopsum[4*pop+3]=kurt;
		}

	return 1;
	}


int Sumstats::getValidSegs(int *nsegs, int npops, int *popsize, int *minind){
	validsegs=ValidSegsVector (nloci);
	for (int locus=0; locus<nloci; ++locus){
		validsegs[locus].assign(nsegs[locus], true);
		for (int seg=0; seg<nsegs[locus]; ++seg){
			int index=seg*npops;
			for (int pop=0; pop<npops; ++pop, ++index){
				if (popsize[pop]-missing[locus][index]<minind[pop]){ validsegs[locus][seg]=false; break; }
				}
			}
		}
	return 1;
	}


int Sumstats::calculateSegs(int *nsegs, int npops, int *popsize, ValidSitesVector &gotValidsites, int scaling){
	validsites=gotValidsites;
	int ncomps=( npops>1 ? (npops*(npops-1))/2 : 1);

	sumstatnames.push_back("Segregating_sites");
	sumstatmap["Segregating_sites"]=OneSumstatVector (nloci*(npops+1), 0);
	sumstatnames.push_back("Segsites_pair");
	sumstatmap["Segsites_pair"]=OneSumstatVector (nloci*ncomps, 0);
	sumstatnames.push_back("Shared_polymorphism");
	sumstatmap["Shared_polymorphism"]=OneSumstatVector (nloci*ncomps, 0);
	sumstatnames.push_back("Fixed_differences");
	sumstatmap["Fixed_differences"]=OneSumstatVector (nloci*ncomps, 0);
	sumstatnames.push_back("Private_pair1");
	sumstatmap["Private_pair1"]=OneSumstatVector (nloci*ncomps, 0);
	sumstatnames.push_back("Private_pair2");
	sumstatmap["Private_pair2"]=OneSumstatVector (nloci*ncomps, 0);
	sumstatnames.push_back("Private_polymorphism");
	sumstatmap["Private_polymorphism"]=OneSumstatVector (nloci*npops, 0);
	sumstatnames.push_back("Theta_Watterson");
	sumstatmap["Theta_Watterson"]=OneSumstatVector (nloci*(npops+1), 0);
	sumstatnames.push_back("Df");
	sumstatmap["Df"]=OneSumstatVector (nloci*ncomps, 0);

	OneSumstatVector &segsites=sumstatmap["Segregating_sites"];
	OneSumstatVector &segpair=sumstatmap["Segsites_pair"];
	OneSumstatVector &shared=sumstatmap["Shared_polymorphism"];
	OneSumstatVector &fixeddiff=sumstatmap["Fixed_differences"];
	OneSumstatVector &privatepair1=sumstatmap["Private_pair1"];
	OneSumstatVector &privatepair2=sumstatmap["Private_pair2"];
	OneSumstatVector &privpol=sumstatmap["Private_polymorphism"];
	OneSumstatVector &theta_w=sumstatmap["Theta_Watterson"];
	OneSumstatVector &df=sumstatmap["Df"];

	int totsize=sum<int>(npops, popsize);
	vector<double> thetawdenom(totsize+1, 0);
	double tempdenom=0;
	for (int i=1; i<totsize; ++i){ tempdenom+=1/(double)i; thetawdenom[i+1]=tempdenom; }

	for (int locus=0; locus<nloci; ++locus){
		for (int seg=0; seg<nsegs[locus]; ++seg){
			if (!validsegs[locus][seg]) continue;	// skip segsite if not covered properly in all populations.
			bool ancestral=false;
			bool derived=false;
			std::vector<short unsigned int> monoancder(npops, 0);	// 0=monomorphic ancestral, 1=polymorphic, 2=monomorphic derived
			int totn=0;
			int index=seg*npops;
			int pop1index=locus*(npops+1);
			for (int pop=0; pop<npops; ++pop, ++index, ++pop1index){
				int validind=popsize[pop]-missing[locus][index]; totn+=validind;
				if (allelemap[locus][index]==0){ monoancder[pop]=0; ancestral=true; }
				else if (allelemap[locus][index]==validind){ monoancder[pop]=2; derived=true; }
				else {
					monoancder[pop]=1;	// both alleles present in population.
					ancestral=true; derived=true;
					++segsites[pop1index];	// number of segregating sites per population.
					theta_w[pop1index]+=1/thetawdenom[validind];
					}
				}
			if (ancestral && derived){
				++segsites[pop1index]; // total number of segregating sites.
				theta_w[pop1index]+=1/thetawdenom[totn];
				}
			for (int pop1=0; pop1<npops; ++pop1){
				if (monoancder[pop1]==1){	// SNP polymorphic in population 1.
					bool polym=false;
					for (int pop2=0; pop2<npops; ++pop2){
						if (pop1 == pop2) continue;
						else if (monoancder[pop2]==1){	// another population is also polymorphic.
							polym=true;
							break;
							}
						}
					if (!polym) ++privpol[locus*npops+pop1]; // number of private polymorphism per population.
					}
				}

			// do pairwise comparisons between populations:
			int comp=0;
			int compindex=locus*ncomps;
			for (std::vector<short unsigned int>::iterator popit1=monoancder.begin(), popend=monoancder.end(); popit1 != popend-1; ++popit1){
				for (std::vector<short unsigned int>::iterator popit2=popit1+1; popit2 != popend; ++popit2, ++comp, ++compindex){
					if (*popit1!=*popit2 || *popit1==1 || *popit2==1){	// site is polymorphic in pair.
						++segpair[compindex]; // number of segregating sites per pair of populations.
						if (*popit1==1 && *popit2==1) ++shared[compindex]; // number of shared polymorphism per pair of populations.
						else if ( (*popit1==0 && *popit2==2) || (*popit1==2 && *popit2==0) ) ++fixeddiff[compindex]; // number of fixed differences per pair of populations.
						else if (*popit1==1 && *popit2!=1) ++privatepair1[compindex]; // number of private polymorphism per pair of groups.
						else if (*popit2==1 && *popit1!=1) ++privatepair2[compindex]; // number of private polymorphism per pair of groups.
						else cout<<"Invalid site pattern in locus: "<<locus<<", segsite: "<<seg<<"!"<<endl;
						}
					} 
				}
			}

		// divide polymorphism counts by number of segregating site per comparison:
		for (int compindex=locus*ncomps; compindex<locus*ncomps+ncomps; ++compindex){
			shared[compindex] /= segpair[compindex];
			df[compindex] = fixeddiff[compindex];
			fixeddiff[compindex] /= segpair[compindex];
			privatepair1[compindex] /= segpair[compindex];
			privatepair2[compindex] /= segpair[compindex];
			}
		for (int pop=0; pop<npops; ++pop) privpol[locus*npops+pop] /= segsites[locus*(npops+1)+npops];

		// divide theta and df estimate by the number of valid sites:
		int nvalidsites=validsites[locus][npops];	// only consider sites valid in all populations.
		if (!nvalidsites){ cout<<"No valid sites for locus "<<locus<<"!"<<endl; nvalidsites=1; }	// what to do with loci without valid sites???
		double scalefactor=(double)scaling/(double)nvalidsites;
		for (int index=locus*(npops+1); index<=locus*(npops+1)+npops; ++index) theta_w[index]*=scalefactor;
		for (int index=locus*ncomps; index<locus*ncomps+ncomps; ++index) df[index]*=scalefactor;
		}

	moments("Segregating_sites", npops+1, segsites.begin(), segsites.end() );
	moments("Segsites_pair", ncomps, segpair.begin(), segpair.end() );
	moments("Shared_polymorphism", ncomps, shared.begin(), shared.end() );
	moments("Fixed_differences", ncomps, fixeddiff.begin(), fixeddiff.end() );
	moments("Private_pair1", ncomps, privatepair1.begin(), privatepair1.end() );
	moments("Private_pair2", ncomps, privatepair2.begin(), privatepair2.end() );
	moments("Private_polymorphism", npops, privpol.begin(), privpol.end() );
	moments("Theta_Watterson", npops+1, theta_w.begin(), theta_w.end() );
	moments("Df", ncomps, df.begin(), df.end() );

	return 1;
	}


int Sumstats::calcDiversity(int *nsegs, int npops, int *popsize, int *minind, ValidSitesVector &gotValidsites, int scaling){
	validsites=gotValidsites;
	int ncomps=( npops>2 ? ((npops*(npops-1))/2)+1 : 1);

	sumstatnames.push_back("Theta_Pi");
	sumstatmap["Theta_Pi"]=OneSumstatVector (nloci*(npops+1), 0);
	sumstatnames.push_back("Shet");
	sumstatmap["Shet"]=OneSumstatVector (nloci*(npops+1), 0);
	sumstatnames.push_back("Dxy");
	sumstatmap["Dxy"]=OneSumstatVector (nloci*ncomps, 0);
	sumstatnames.push_back("Fst");
	sumstatmap["Fst"]=OneSumstatVector (nloci, 0);
	sumstatnames.push_back("Phist");
	sumstatmap["Phist"]=OneSumstatVector (nloci*ncomps, 0);

	OneSumstatVector &segsites=sumstatmap["Segregating_sites"];
	OneSumstatVector &thetap=sumstatmap["Theta_Pi"];
	OneSumstatVector &shet=sumstatmap["Shet"];
	OneSumstatVector &pdxy=sumstatmap["Dxy"];
	OneSumstatVector &Fst=sumstatmap["Fst"];
	OneSumstatVector &Phist=sumstatmap["Phist"];

	for (int locus=0; locus<nloci; ++locus){
		vector<double> pi(npops+1, 0);
		vector<double> dxy(ncomps, 0);
		vector<double> sumvart(ncomps, 0);
		vector<double> sumvara(ncomps, 0);
		for (int seg=0; seg<nsegs[locus]; ++seg){
			if (validsegs[locus][seg]){	// skip segsite if not covered properly in all populations.
				double ssdwp=0;
				int totder=0, totanc=0, totn=0, pair=0, valid=1;
				for (int pop1=0, index1=seg*npops; pop1<npops; ++pop1, ++index1){
					int validind1=popsize[pop1]-missing[locus][index1];
					int denom1=validind1*(validind1-1);
					totn+=validind1;
					int nderived1=allelemap[locus][index1];
					totder+=nderived1;
					int nancestral1=validind1-nderived1;
					totanc+=nancestral1;
//					int dersquared1=nderived1*(nderived1-1);
//					int ancsquared1=nancestral1*(nancestral1-1);
//					double total1=(dersquared1+ancsquared1)/(double)denom1;
					double diffprop1=(2*nderived1*nancestral1)/(double)denom1;
					ssdwp+=((double)(validind1-1)/2)*diffprop1;
					pi[pop1]+=diffprop1;
					for (int pop2=pop1+1, index2=seg*npops+pop1+1; pop2<npops; ++pop2, ++pair, ++index2){
						int validind2=popsize[pop2]-missing[locus][index2];
						int validpair=validind1+validind2;
						int denom2=validind2*(validind2-1);
						int denompair=validind1*validind2;
						int denomtot=validpair*(validpair-1);
						int nderived2=allelemap[locus][index2];
						int nancestral2=validind2-nderived2;
//						int dersquared2=nderived2*(nderived2-1);
//						int dersquaredpair=nderived1*nderived2;
//						int dersquaredtot=(nderived1+nderived2)*(nderived1+nderived2-1);
//						int ancsquared2=nancestral2*(nancestral2-1);
//						int ancsquaredpair=nancestral1*nancestral2;
//						int ancsquaredtot=(nancestral1+nancestral2)*(nancestral1+nancestral2-1);
//						double total2=(dersquared2+ancsquared2)/(double)denom2;
//						double totalpair=(dersquaredpair+ancsquaredpair)/(double)denompair;
//						double totaltot=(dersquaredtot+ancsquaredtot)/(double)denomtot;
						double diffprop2=(2*nderived2*nancestral2)/(double)denom2;
						double diffproppair=(nderived1*nancestral2+nancestral1*nderived2)/(double)denompair;
						double diffproptot=(2*(nderived1+nderived2)*(nancestral1+nancestral2))/(double)denomtot;
						dxy[pair]+=diffproppair;
						double ssdwppair=((double)(validind1-1)/2)*diffprop1+((double)(validind2-1)/2)*diffprop2;
						double ssdtotpair=((double)(validpair-1)/2)*diffproptot;
						double sumnpair=(double)validpair-((validind1*validind1)/(double)validpair+(validind2*validind2)/(double)validpair);
						sumvart[pair]+=(validpair-2)*ssdtotpair-(validpair-1-sumnpair)*ssdwppair;
						sumvara[pair]+=(validpair-2)*ssdtotpair-(validpair-1)*ssdwppair;
						}
					}

//				int totdersquared=totder*(totder-1);
//				int totancsquared=totanc*(totanc-1);
				int denom=totn*(totn-1);
//				double total=(totdersquared+totancsquared)/(double)denom;
				double diffprop=(2*totder*totanc)/(double)denom;
				double ssdtot=((double)(totn-1)/2)*diffprop;
				double sumn=0;
				for (int pop=0; pop<npops; ++pop) sumn+=(pow((double)(popsize[pop]-missing[locus][seg*npops+pop]), (double)2))/(double)totn;
				sumvart[pair]+=(totn-npops)*ssdtot-(totn-1-sumn)*ssdwp;
				sumvara[pair]+=(totn-npops)*ssdtot-(totn-1)*ssdwp;
				pi[npops]+=diffprop;
				}
			}

		double pi_av=0;
		int popindex=locus*(npops+1);
		int nvalidsites=validsites[locus][npops];	// only consider sites valid in all populations.
		if (!nvalidsites){ cout<<"No valid sites for locus "<<locus<<"!"<<endl; nvalidsites=1; }	// what to do with loci without valid sites???
		double scalefactor=(double)scaling/(double)nvalidsites;
		for (int pop=0; pop<npops; ++pop, ++popindex){
			if ( segsites[popindex] ) shet[popindex]=pi[pop] / segsites[popindex];	// calculate site heterozygosity before scaling.
			else shet[popindex]=0;
			pi[pop]*=scalefactor;
			thetap[popindex]=pi[pop];
			pi_av+=pi[pop];
			}
		if (segsites[locus*(npops+1)+npops]) shet[popindex]=pi[npops] / segsites[locus*(npops+1)+npops];
		else shet[popindex]=0;
		pi[npops]*=scalefactor;
		thetap[popindex]=pi[npops];
		pi_av=pi_av/(double)npops;
		Fst[locus]=1-(pi_av/pi[npops]); // only overall Fst implemented! Maybe need to use weighted average of within population pi?
		double dxy_av=0;
		int iter=( npops>2 ? ((npops*(npops-1))/2) : 1);
		int compindex=locus*ncomps;
		for (int i=0; i<iter; ++i, ++compindex){
			dxy[i]*=scalefactor;
			pdxy[compindex]=dxy[i];
			dxy_av+=dxy[i];
			Phist[compindex]=sumvara[i]/sumvart[i];
			}
		if (npops>2){
			pdxy[compindex]=dxy_av/(double)iter;
			Phist[compindex]=sumvara[ncomps-1]/sumvart[ncomps-1];
			}
		}

	moments("Theta_Pi", npops+1, thetap.begin(), thetap.end() );
	moments("Shet", npops+1, shet.begin(), shet.end() );
	moments("Dxy", ncomps, pdxy.begin(), pdxy.end() );
	moments("Fst", 1, Fst.begin(), Fst.end() );
	moments("Phist", ncomps, Phist.begin(), Phist.end() );

	return 1;
	}


double Sumstats::denominator(int n, int segsites){
	double epsilon=1e-10;
	double a1=0, a2=0, b1=0, b2=0, c1=0, c2=0, e1=0, e2=0, denom=0;

	for (int ind=1; ind<n; ind++){
		a1+=(double)1/(double)ind;
		a2+=(double)1/(double)(ind*ind);
		}

	b1=((double)(n+1))/((double)(3*(n-1)));
	b2=((double)(2*(n*n+n+3)))/((double)(9*n*(n-1)));
	c1=b1-((double)1)/a1;

	if(fabs(c1) < epsilon) c1=0;
 	c2=b2-(n+2)/(a1*n)+(a2/(a1*a1));
	if(fabs(c2) < epsilon) c2=0;
	e1=c1/a1;
	if(fabs(e1) < epsilon) e1=0;
	e2=c2/(a1*a1 + a2);
	if(fabs(e2) < epsilon) e2=0;
    
	denom=(e1*segsites + e2*segsites*((double)segsites-1) );
	if(fabs(denom)<epsilon) denom=epsilon;
	if(denom<=0) return 0;
	double denom2=sqrt(denom);
	return denom2;
	}


int Sumstats::tajD(int npops, int *popsize, bool teststat){
	sumstatnames.push_back("Tajimas_D");
	sumstatmap["Tajimas_D"]=OneSumstatVector (nloci*(npops+1), 0);

	OneSumstatVector &segsites=sumstatmap["Segregating_sites"];
	OneSumstatVector &pthetap=sumstatmap["Theta_Pi"];
	OneSumstatVector &pthetaw=sumstatmap["Theta_Watterson"];
	OneSumstatVector &ptajd=sumstatmap["Tajimas_D"];

	int totsize=0;
	for (int locus=0; locus<nloci; locus++){
		for (int pop=0; pop<=npops; pop++){
			int index=locus*(npops+1)+pop;
			int nind=0;
			if (pop<npops){nind=popsize[pop]; totsize+=nind;}
			else nind=totsize;
			int nsegs=segsites[index];
			double thetap=pthetap[index];
			double thetaw=pthetaw[index];
			double tajd;
			if (teststat){
				double denom=denominator(nind, nsegs);
				if (denom<=0) return 0;
				else tajd=(thetap-thetaw)/denom;
				}
			else tajd=thetap-thetaw;
			ptajd[index]=tajd;
			}
		}
	moments("Tajimas_D", npops+1, ptajd.begin(), ptajd.end() );
	return 1;
	}


int Sumstats::ZnS(int *nsegs, int npops, int *popsize, int *minind, char ***haplotypes){
	sumstatnames.push_back("ZnS");
	sumstatmap["ZnS"]=OneSumstatVector (nloci*(npops+1), 0);
	OneSumstatVector &ZnS=sumstatmap["ZnS"];

	for (int pop=0; pop<npops; pop++){
		if (popsize[pop] % 2){ cout<<"Uneven number of individuals per population, assuming haplotized data."<<endl; return 0; }
		}

	int totind=sum(npops, popsize)/2;
	vector<pair<int, int> > geno; geno.reserve(totind);
	double rsqr=0;
//	vector<int> validind(npops, 0);

	for (int locus=0; locus<nloci; ++locus){
		int segpairs[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
		for (int seg1=0; seg1<nsegs[locus]-1; ++seg1){
			for (int seg2=seg1+1; seg2<nsegs[locus]; ++seg2){
				int allind=0, startindex=0;
				bool validpair=true;
				for (int pop=0, index=locus*(npops+1); pop<npops; ++pop, ++index){
					int validind=0;
					for (int ind=0; ind<popsize[pop]; ind+=2){
						char *indhap1=haplotypes[locus][2*allind];
						char *indhap2=haplotypes[locus][2*allind+1];
						++allind;
						if (indhap1[seg1]=='2' || indhap2[seg1]=='2') continue;
						else if (indhap1[seg2]=='2' || indhap2[seg2]=='2') continue;
						else {
							geno.push_back(pair<int, int>(indhap1[seg1]+indhap2[seg1], indhap1[seg2]+indhap2[seg2]));
							++validind;
							}
						}
					if (validind>=minind[pop]/2){ rsqr=r2( geno.begin()+startindex, geno.end() ); if ( isfinite(rsqr) ){ ZnS[index]+=rsqr; ++segpairs[pop]; }; }
					else validpair=false;
					startindex+=validind;
					}
				if (validpair){ rsqr=r2( geno.begin(), geno.end() ); if ( isfinite(rsqr) ){ ZnS[locus*(npops+1)+npops]+=rsqr; ++segpairs[npops]; }; }
				geno.clear();
				}
			}
		for (int pop=0, index=locus*(npops+1); pop<=npops; ++pop, ++index){ if (segpairs[pop]) ZnS[index] /= segpairs[pop]; else ZnS[index]=nan(""); }
		}
	moments("ZnS", npops+1, ZnS.begin(), ZnS.end() );
	return 1;
	}


int Sumstats::ZnS_Quartiles(int *nsegs, int npops, int *popsize, int *minind, char ***haplotypes){
	sumstatnames.push_back("ZnS_quart");
	sumstatmap["ZnS_quart"]=OneSumstatVector (nloci*(npops+1), 0);
	OneSumstatVector &ZnS=sumstatmap["ZnS_quart"];

	for (int pop=0; pop<npops; pop++){
		if (popsize[pop] % 2){ cout<<"Uneven number of individuals per population, assuming haplotized data."<<endl; return 0; }
		}

	int totind=sum<int>(npops, popsize)/2;
	vector<pair<int, int> > geno; geno.reserve(totind);
	double rsqr=0;
//	vector<int> validind(npops, 0);

	for (int locus=0; locus<nloci; locus++){
		int segpairs[]={0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
//		cout<<"Locus: "<<locus<<" - "<<nsegs[locus]<<"/"<<quartiles[3*locus]<<"/"<<quartiles[3*locus+2]<<endl;
		for (int seg1=0; seg1<quartiles[3*locus]; ++seg1){
			for (int seg2=quartiles[3*locus+2]; seg2<nsegs[locus]; ++seg2){
				int allind=0, startindex=0;
				bool validpair=true;
				for (int pop=0; pop<npops; pop++){
					int index=locus*(npops+1)+pop;
					int validind=0;
					for (int ind=0; ind<popsize[pop]; ind+=2){
						char *indhap1=haplotypes[locus][2*allind];
						char *indhap2=haplotypes[locus][2*allind+1];
						allind++;
						if (indhap1[seg1]=='2' || indhap2[seg1]=='2') continue;
						else if (indhap1[seg2]=='2' || indhap2[seg2]=='2') continue;
						else {
							geno.push_back(pair<int, int>(indhap1[seg1]+indhap2[seg1], indhap1[seg2]+indhap2[seg2]));
							validind++;
							}
						}
					if (validind>=minind[pop]/2){ rsqr=r2( geno.begin()+startindex, geno.end() ); if ( isfinite(rsqr) ){ ZnS[index]+=rsqr; segpairs[pop]++; }; }
					else validpair=false;
					startindex+=validind;
					}
				if (validpair){ rsqr=r2( geno.begin(), geno.end() ); if ( isfinite(rsqr) ){ ZnS[locus*(npops+1)+npops]+=rsqr; segpairs[npops]++; }; }
				geno.clear();
				}
			}
		for (int pop=0; pop<=npops; pop++) ZnS[locus*(npops+1)+pop] /= segpairs[pop];
		}
	moments("ZnS_quart", npops+1, ZnS.begin(), ZnS.end() );
	return 1;
	}


inline double Sumstats::r2(vector<pair<int, int> >::iterator start, vector<pair<int, int> >::iterator end){
	double avgx=0, avgy=0, sumxx=0, sumyy=0, sumxy=0;
	int n=0;

	for (vector<pair<int, int> >::iterator it=start, ite=end; it != ite; ++it){
		avgx+=it->first;
		avgy+=it->second;
		++n;
		}

	avgx /= n;
	avgy /= n;

	for (vector<pair<int, int> >::iterator it=start, ite=end; it != ite; ++it){
		double dx=it->first-avgx;
		double dy=it->second-avgy;
		sumxx+=dx*dx;
		sumyy+=dy*dy;
		sumxy+=dx*dy;
		}

	return pow(sumxy/sqrt(sumxx*sumyy), 2);
	}


int Sumstats::sfs(bool folded, int *nsegs, int npops, int *popsize, int *minind){
	if ( npops > 2 ){ cout<<"Site frequency spectrum cannot be calculated for more than two populations!\n"; return 0; }
	int ncomps=0, nder[2], halfsize[2];
	for (int pop=0; pop<npops; ++pop) halfsize[pop]=popsize[pop]/2+1;
	if (folded) ncomps=( npops==2 ? halfsize[0]*halfsize[1] : halfsize[0]);
	else ncomps=( npops==2 ? (popsize[0]+1)*(popsize[1]+1) : popsize[0]+1);

	sumstatnames.push_back("SFS");
	sumstatmap["SFS"]=OneSumstatVector (ncomps, 0);
	OneSumstatVector &psfs=sumstatmap["SFS"];

	for (int locus=0; locus<nloci; ++locus){
		for (int seg=0; seg<nsegs[locus]; ++seg){
			for (int pop=0; pop<npops; ++pop){
				int index=seg*npops+pop;
				int nmissing=missing[locus][index];
				int validind=popsize[pop]-nmissing;
				if (validind>=minind[pop]){
					nder[pop]=allelemap[locus][index];
					if (nmissing) for (int draw=0; draw<nmissing; ++draw){
						double fder=allelemap[locus][index] / (double)validind;
						if ( ((double)rand() / (double)RAND_MAX) < fder ) ++nder[pop];
						}
					if (folded && nder[pop] >= halfsize[pop] ) nder[pop]=popsize[pop]-nder[pop];
					}
				}
			if (npops==1) ++psfs[ nder[0] ];
			else ++psfs[ nder[0]*halfsize[1]+nder[1] ];
			}
		}

	locsummarymap["SFS"]=std::vector<double>(ncomps, 0);
	std::vector<double> &itpop=locsummarymap["SFS"];
	for (int comp=0; comp<ncomps; ++comp) itpop[comp]=psfs[comp] / nloci;

	return 1;
	}


int Sumstats::printSumstats(char *outfile){
	ofstream output (outfile);
	if ( !output.is_open() ){
		cout<<"File "<<outfile<<" could not be opened!\n";
		return 0;
		}

	int nsum=sumstatnames.size();
	cout<<"nsum: "<<nsum<<endl;
	for (int sumstat=0; sumstat<nsum; ++sumstat){
		int npops=locsummary[sumstat].size();
		if (sumstatnames[sumstat] != "SFS"){
			for (int pop=0; pop<npops/2; ++pop) output<<sumstatnames[sumstat]<<"_"<<pop<<"_mean\t"<<
				sumstatnames[sumstat]<<"_"<<pop<<"_sd\t";
			}
		else for (int comp=0; comp<npops; ++comp) output<<sumstatnames[sumstat]<<"_"<<comp<<"\t";
		}
	output<<endl;
	for (int sumstat=0; sumstat<nsum; ++sumstat){
		int npops=locsummary[sumstat].size();
		if (sumstatnames[sumstat] != "SFS"){
			for (int pop=0; pop<npops; pop+=2) output<<locsummary[sumstat][pop]<<"\t"<<locsummary[sumstat][pop+1]<<"\t";
			}
		else for (int comp=0; comp<npops; ++comp) output<<locsummary[sumstat][comp]<<"\t";
		}
	output<<endl;
	output.close();
	return 1;
	}

int Sumstats::printSumstatsPriors(char *outfile, int nsim, TPriorVector *priors, bool printHeader, bool append){
	ofstream outstream;
	if (append) outstream.open(outfile, ios::out | ios::app);
	else outstream.open(outfile);
	if ( !outstream.is_open() ){
		cout<<"File "<<outfile<<" could not be opened!\n";
		return 0;
		}

	int nsimplepriors=(priors!=NULL) ? priors->numSimplePrior : 0;
	int ncombpriors=(priors!=NULL) ? priors->numCombinedPrior : 0;
	int nsum=sumstatnames.size();
//	cout<<"nsum: "<<nsum<<endl;
	if (printHeader){
		outstream<<"NSim\t";
		for (int prior=0; prior<nsimplepriors; ++prior) if (priors->simplePriors[prior]->output) outstream<<priors->simplePriors[prior]->name<<"\t";
		for (int prior=0; prior<ncombpriors; ++prior) if (priors->combinedPriors[prior]->output) outstream<<priors->combinedPriors[prior]->name<<"\t";
		for (std::vector<std::string>::iterator ssname=sumstatnames.begin(); ssname!=sumstatnames.end(); ++ssname){
			int npops=locsummarymap[*ssname].size();
			if (*ssname != "SFS"){
				for (int pop=0; pop<npops/4; ++pop) outstream<<*ssname<<"_"<<pop<<"_mean\t"<<
					*ssname<<"_"<<pop<<"_var\t"<<
					*ssname<<"_"<<pop<<"_skew\t"<<
					*ssname<<"_"<<pop<<"_kurt\t";
				}
			}
		outstream<<endl;
		}
	outstream<<nsim<<"\t";
	for (int prior=0; prior<nsimplepriors; ++prior) if (priors->simplePriors[prior]->output) outstream<<priors->simplePriors[prior]->curValue<<"\t";
	for (int prior=0; prior<ncombpriors; ++prior) if (priors->combinedPriors[prior]->output) outstream<<priors->combinedPriors[prior]->curValue<<"\t";
	for (std::vector<std::string>::iterator ssname=sumstatnames.begin(); ssname!=sumstatnames.end(); ++ssname){
		int npops=locsummarymap[*ssname].size();
		if (*ssname != "SFS"){
			for (int pop=0; pop<npops; pop+=4) outstream<<locsummarymap[*ssname][pop]<<"\t"<<locsummarymap[*ssname][pop+1]<<"\t"<<locsummarymap[*ssname][pop+2]<<"\t"<<locsummarymap[*ssname][pop+3]<<"\t";
			}
		}
	outstream<<endl;
	outstream.close();
	return 1;
	}





