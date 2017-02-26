#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <time.h>
#include "hyperpriors.h"
#include "mssim.h"
#include "input.h"
#include "missing.h"
#include "sumstats.h"
#include "TException.h"
#include "TLog.h"
#include "TRandomGenerator.h"
#include "TPrior.h"
#include "misc.h"

#define VERSION "0.50"

using namespace std;



int main(int argc, char **argv) {
	cout<<"postABC version "<<VERSION<<endl;
	clock_t starttime=clock();
	char *infile=NULL, *outfile=NULL, *misslist=NULL, *validsites=NULL;
	string estfile;
	char **params, **msargs;
	int seed=0, numsim=1, npars=0, nmsargs=0, nloci=1, nref=10000, length=2000, scaling=1, npops=0, submode=0;
	int popsize[10];
	int minind[]={1,1,1,1,1,1,1,1,1,1};
	int subind[10];
	bool usesupermatrix=false, printsfs=false, folded=true, subsamp=false, ssbefore=true, haplotize=false, persitecalc=false, calcZnS=false, calcZnS_quart=false;

	if (argc<9) {
		cout<<"Usage is <npops> --seed <seed for random number generator> --infile <filename of output file from ms> --output <filename for summary statistics output> --missing <filename of missing file list> --hypr <number of hyperpriors> [ <type of distribution> <mean> <shape> <left limit> <right limit> for each hyperprior ] --ms <ms command line> --popsizes <nind pop1> <nind pop2> ... --minind <minimal number of individuals for sumstat calclation pop1> <pop2> ... --nloci <number of loci>\n";
		exit(1);
	} else {
		npops=atoi(argv[1]);
		for (int i=2; i<argc; ++i) {
			if (i+1 != argc){
				if (string(argv[i]) == "--seed") seed=atoi(argv[++i]);
				else if (string(argv[i]) == "--infile") infile=argv[++i];
				else if (string(argv[i]) == "--estfile") estfile=string(argv[++i]);
				else if (string(argv[i]) == "--output") outfile=argv[++i];
				else if (string(argv[i]) == "--numsim") numsim=atoi(argv[++i]);
				else if (string(argv[i]) == "--hypr"){
					npars=atoi(argv[++i]);
					params=argv+i+1; i+=(5*npars); }
				else if (string(argv[i]) == "--ms"){
					msargs=argv+i; ++nmsargs;
					while ( strlen(argv[i+1]) == 1 || ( strlen(argv[i+1]) > 1 && argv[i+1][1] != '-' ) ){ ++nmsargs; ++i; }
					}
				else if (string(argv[i]) == "--missing") misslist=argv[++i];
				else if (string(argv[i]) == "--popsizes") for (int j=0; j<npops; ++j) popsize[j]=atoi(argv[++i]);
				else if (string(argv[i]) == "--nloci") nloci=atoi(argv[++i]);
				else if (string(argv[i]) == "--loclength") length=atoi(argv[++i]);
				else if (string(argv[i]) == "--scaling") scaling=atoi(argv[++i]);
				else if (string(argv[i]) == "--validsites"){
					validsites=argv[++i];
					persitecalc=true;
					}
				else if (string(argv[i]) == "--Nref") nref=atoi(argv[++i]);
				else if (string(argv[i]) == "--minind") for (int j=0; j<npops; ++j) minind[j]=atoi(argv[++i]);
				else if (string(argv[i]) == "--subsample"){
					subsamp=true;
					submode==atoi(argv[++i]);
					haplotize=( 1 == atoi(argv[++i]) );
					for (int j=0; j<npops; ++j) subind[j]=atoi(argv[++i]); }
				else if (string(argv[i]) == "--supermatrix") usesupermatrix=( 1 == atoi(argv[++i]) );
				else if (string(argv[i]) == "--sfs") printsfs=( 1 == atoi(argv[++i]) );
				else if (string(argv[i]) == "--ZnS") calcZnS=( 1 == atoi(argv[++i]) );
				else if (string(argv[i]) == "--ZnS_quart") calcZnS_quart=( 1 == atoi(argv[++i]) );
				else cout<<"Unkown parameter "<<argv[i]<<endl;
				}
			}
		}

	if (npops>10){
		cout<<"More than 10 populations are currently not supported!"<<endl;
		exit(1);
		}
	int totsize=sum<int>(npops, popsize);
	if (nmsargs){
		if ( totsize!=atoi(msargs[1]) ){ cout<<"Wrong number of samples specified in MS command line!"<<endl; exit(1); }
		if ( nloci!=atoi(msargs[2]) ){ cout<<"Number of loci specified does not correspond to number of loci specified in ms parameters!"<<endl; exit(1); }
		}
	if (!subsamp) for (int index=0; index<npops; ++index) subind[index]=popsize[index];
	seed==0 ? srand(time(NULL)) : srand(seed);
	clock_t opttime=clock();
	cout<<"Options successfully processed."<<endl;

//******************************************************************************************************************

	// initialize random number generator:
	TLog *logFile=new TLog;
	TRandomGenerator *randomGenerator=new TRandomGenerator(seed, true);	// initialize with fixed seed.
	cout<<"Random number generator successfully initialized."<<endl;

	// read file with prior distributions:
	TPriorVector *priors=NULL;
	if ( !estfile.empty() ){
		try {
			priors=new TPriorVector(estfile, logFile, randomGenerator);
			}
		catch (TException & error){
			cout<< error.getMessage() <<endl;
			}
		}
	char **SimValues=allocation<char>(nmsargs, 30);
//	for(int i=0; i<nmsargs; ++i) strcpy(SimValues[i], msargs[i]);	// make copy of msarg array of c-strings.
	char **HypValues=allocation<char>(npars*5, 30);
	for(int i=0; i<(npars*5); ++i) strcpy(HypValues[i], params[i]);	// make copy of hypparam array of c-strings.

	// read files with missing data:
	MissingData missdata(nloci, totsize);
	if ( misslist != NULL ){
		if ( !missdata.ReadMisslist(misslist) ) exit(1);
		else cout<<"Processing missing data..."<<endl;
		if ( !missdata.ReadMissing() ) exit(1);
		if (usesupermatrix) missdata.createMissMatrix(length);
		}
	if ( !missdata.getValidSites(length, npops, popsize, minind) ) exit(1);

	cout<<"Starting the simulation loops..."<<endl;

	// loop through all simulation iterations:
	for (int sim=0; sim<numsim; ++sim){
		clock_t loopstarttime=clock();
		if (priors!=NULL) priors->getNewValues();	// create new random values for each parameter.


		if( !estfile.empty() ){
			// parse arrays of character pointers for simulation parameters:
			for(int i=0; i<nmsargs; ++i){
				strcpy(SimValues[i], msargs[i]);		// ms changes the tbs tags with actual values from the hyperprior distribution!
				if ( priors->isPrior( string(msargs[i]) ) ){
//					SimValues[i][0]='\0';	// set c-string to empty.
					float val=priors->getValue( string(msargs[i]) );
					cout<<"Identified prior tag "<<msargs[i]<<". Replaced tag with value: "<<val<<"."<<endl;
					if (val!=_nan) snprintf(SimValues[i], 30, "%f", val);
					}
				}
			// parse arrays of character pointers for hyperparameters:
			for(int i=0; i<(npars*5); ++i){
				if ( priors->isPrior( string(params[i]) ) ){
					HypValues[i][0]='\0';	// set c-string to empty.
					float val=priors->getValue( string(params[i]) );
					cout<<"Identified hyperprior tag "<<params[i]<<". Replaced tag with value: "<<val<<"."<<endl;
					if (val!=_nan) snprintf(HypValues[i], 30, "%f", val);
					}
				}
			}

		cout<<"Simulation "<<sim+1<<" of "<<numsim<<"."<<endl;
		cout<<endl<<endl;
		cout<<"Simulation parameters:"<<endl;
		for(int i=0; i<nmsargs; ++i) cout<<SimValues[i]<<" ";
		cout<<endl<<endl;
		cout<<"Hyperprior parameters:"<<endl;
		for(int i=0; i<(npars*5); ++i) cout<<HypValues[i]<<" ";
		cout<<endl<<endl;

//******************************************************************************************************************

		RawData simdata(nloci, npops, popsize);

		if (nmsargs){
			HyperPriors hyperpr(npars, nloci);
			hyperpr.hyperpriors(nref, length, HypValues);
			if ( mssim(nmsargs, SimValues, npars, hyperpr.dist, simdata.nsegs, simdata.positions, simdata.haplotypes)==nloci ) cout<<"Data successfully simulated."<<endl;
			else { cout<<"Error in MS data simulation step."<<endl; exit(1); }
			}
		else {
			int lociinfile=simdata.readms_c(infile);
			if (!lociinfile){ cout<<"Could not acquire locus data from ms file!"<<endl; exit(1); }
			else if (lociinfile!=nloci){ cout<<"Number of loci specified in command line does not match number of loci in input file!"<<endl; exit(1); }
			else cout<<"MS file successfully read in."<<endl;
			}
		clock_t readtime=clock();

		Sumstats sumstats(nloci);
		if ( simdata.processms(missdata, nloci, subsamp, ssbefore, submode, haplotize, subind, sumstats) ) cout<<"Simulation data successfully processed."<<endl;
		else { cout<<"Could not process simulation data!"<<endl; exit(1); }
		if ( sumstats.getValidSegs(simdata.nsegs, npops, popsize, minind) ) cout<<"Valid segregating sites successfully processed."<<endl;
		else { cout<<"Could not determine valid segregating sites!"<<endl; exit(1); }
		clock_t processtime=clock();

		if (persitecalc){
			if ( sumstats.readgroupinfo(validsites) ) cout<<"File with valid sites information successfully processed."<<endl;
			else { cout<<"Could not read file with valid sites information!"<<endl; exit(1); }
			}
		clock_t nvalidtime=clock();

		sumstats.calculateSegs(simdata.nsegs, npops, subind, missdata.validsites, scaling);
		sumstats.calcDiversity(simdata.nsegs, npops, subind, minind, missdata.validsites, scaling);
		sumstats.tajD(npops, subind, false);
		if (calcZnS || calcZnS_quart){
			if (!haplotize && submode==0){
				if (calcZnS) sumstats.ZnS(simdata.nsegs, npops, subind, minind, simdata.subhaplotypes);
				else if (calcZnS_quart) sumstats.ZnS_Quartiles(simdata.nsegs, npops, subind, minind, simdata.subhaplotypes);
				}
			else cout<<"Cannot calculate ZnS statistic for haplotized data!"<<endl;
			}
		if (printsfs) sumstats.sfs(folded, simdata.nsegs, npops, subind, minind);
		clock_t sumstattime=clock();
		cout<<"Sumstats successfully calculated."<<endl;

		if (!sim) sumstats.printSumstatsPriors(outfile, sim, priors, 1, 0);	// first simulation round, print header with output.
		else sumstats.printSumstatsPriors(outfile, sim, priors, 0, 1);
		clock_t printtime=clock();
		cout<<"Sumstats successfully printed."<<endl;
		clock_t loopendtime=clock();
		printf ("Finished loop %d, execution time: %ld clicks (%f seconds).\n", loopendtime-loopstarttime,((float)(loopendtime-loopstarttime))/CLOCKS_PER_SEC);
		}

	// free memory:
	freemem<char>(SimValues, nmsargs);
	freemem<char>(HypValues, npars*5);
	delete logFile;
	delete randomGenerator;
	delete priors;

	clock_t endtime=clock();
//	printf ("Option processing time: %ld clicks (%f seconds).\n", opttime-starttime,((float)(opttime-starttime))/CLOCKS_PER_SEC);
//	printf ("Input reading time: %ld clicks (%f seconds).\n", readtime-opttime,((float)(readtime-opttime))/CLOCKS_PER_SEC);
//	printf ("Data processing time: %ld clicks (%f seconds).\n", processtime-readtime,((float)(processtime-readtime))/CLOCKS_PER_SEC);
//	printf ("Validsites reading time: %ld clicks (%f seconds).\n", nvalidtime-processtime,((float)(nvalidtime-processtime))/CLOCKS_PER_SEC);
//	printf ("Summary statistics calculation time: %ld clicks (%f seconds).\n", sumstattime-nvalidtime,((float)(sumstattime-nvalidtime))/CLOCKS_PER_SEC);
//	printf ("Memory freeing time: %ld clicks (%f seconds).\n", endtime-printtime,((float)(endtime-printtime))/CLOCKS_PER_SEC);
	printf ("Execution time: %ld clicks (%f seconds).\n", endtime,((float)endtime)/CLOCKS_PER_SEC);

	return 0;
	}


