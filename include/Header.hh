#ifndef GeneralHeader
#define GeneralHeader

using namespace std;
//using namespace __gnu_cxx;

#include <stdio.h>
#include <stdlib.h>
#include <ext/numeric>
#include <cmath>
#include <list>
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <map>
#include <iterator>
#include <algorithm>
#include <boost/utility.hpp>
#include <cstring>
#include <sstream>
#include <fstream>
#include "dSFMT.h"	// Kirsten used it for some functions, so lets just stick with it
//#include "Xmgrace.hh"	// Grace pipeline added by Brem
#include <sys/time.h>
#include <ctime>
#include <typeinfo>
#include <string>

#define toDigit(c) (c-'0')  // Converts char to digit

//Constants for Gene.hh
const int binding_length = 10;

//Constants defined here, used in Population.hh
//Currently gives a Segmentation fault if NR and NC are not equal!
const int NR=100;  //While there is no diffusion, there seems to be no reason to make the field too big.
const int NC=100;

//Constants defined here, used in Genome.cc
const int init_nr_gene_types = 10;
const int init_nr_tfbs_per_gene = 3;
const int WeightRange = 3;  //TFBS weights can range from -WeightRange to +WeightRange.

// Mutations
const double gene_threshold_mu = 0.001;
const double gene_activity_mu = 0.001;
const double gene_binding_domain_mu = 0.001;

const double tfbs_binding_site_mu = 0.001;
const double tfbs_activity_mu = 0.001;

const double gene_duplication_mu = 0.001;
const double gene_deletion_mu = 0.001;

const double tfbs_duplication_mu = 0.001;
const double tfbs_deletion_mu = 0.001;

const int tfbs_selection_exponent = 10;
const double empty_tf_claim_zero = 1.0;

//constants for Population.cc
const int SimTime=1000;
const int TimeSaveGrid = 100; //How many timesteps to save the whole grid as raw data.
const int TimeTerminalOutput = 10;  //Note that this is also the check for extinction, so let it check!
const double death_rate = 0.03;
const double repl_rate = 1.0;
const int replication_neighbourhood = 3;  //i.e. a 3x3 grid represents the neighbourhood.

//constants for Prokaryote.cc

//Variables defined in World.cc
extern dsfmt_t dsfmt;
inline double uniform() { return dsfmt_genrand_close_open(&dsfmt); }

extern int Time;
extern int initial_seed;
extern int TargetExpression[5];
extern string folder;
extern bool mutational_neighbourhood;
extern int NrMutants;
extern string genome_init;
extern string genestate_init;
extern bool mutations_on;

const string genome_file="/home/sam/Documents/Endosymbiosis/Model/Prokaryotes1.0/input/Caulobacter_SO_genome.g";
//const string genome_file="/home/sam/Documents/Endosymbiosis/Model/Projects/Caulobacter_Alpha/S107/MA_100k.g";
const string genestate_file="/home/sam/Documents/Endosymbiosis/Model/Prokaryotes1.0/input/M_stage_genestates.g";
//const string genestate_file="/home/sam/Documents/Endosymbiosis/Model/Projects/Caulobacter_Alpha/S107/MA_100k_MA_expr.g";
//const string genome_file="";
//const string genestate_file="";


const bool StageTargets[4][5] = {
  false, false, true, false, true,      // 0 0 1 0 1    G1
  false, true, false, false, false,     // 0 1 0 0 0    S
  true, false, false, false, false,     // 1 0 0 0 0    G2
  true, false, false, true, true        // 1 0 0 1 1    M
};

#endif
