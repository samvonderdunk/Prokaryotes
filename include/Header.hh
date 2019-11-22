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
const int repl_step_size = 20;

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
const int TimeTerminalOutput = 1;  //Note that this is also the check for extinction, so let it check!
const int TimeSaveGrid = 100; //How many timesteps to save the whole grid as raw data.
const int TimeSaveBackup = 100;
const int TimePruneFossils = 100;
const int TimeOutputFossils = 100;
const double death_rate = 0.02;
const int TimeZero=0;
const double repl_rate = 1.0;
const int replication_neighbourhood = 3;  //i.e. a 3x3 grid represents the neighbourhood.
const int generation_sample = 10000;

//constants for Prokaryote.cc

//Variables defined in World.cc
extern dsfmt_t dsfmt;
inline double uniform() { return dsfmt_genrand_close_open(&dsfmt); }

extern int Time;
extern int initial_seed;
extern string folder;
extern bool mutational_neighbourhood;
extern bool attractor_landscape;
extern int NrMutants;
extern int NrInitialStates;
extern string genome_init;
extern string genestate_init;
extern string backup_reboot;
extern string anctrace_reboot;
extern bool mutations_on;

const string genome_file="/home/sam/Documents/Endosymbiosis/Model/Prokaryotes1.0/input/Caulobacter_SO_genome.g";
//const string genome_file="/home/sam/Documents/Endosymbiosis/Model/Projects/Caulobacter_Alpha/S107/MA_100k.g";
const string genestate_file="/home/sam/Documents/Endosymbiosis/Model/Prokaryotes1.0/input/M_stage_genestates.g";
//const string genestate_file="/home/sam/Documents/Endosymbiosis/Model/Projects/Caulobacter_Alpha/S107/MA_100k_MA_expr.g";
// const string backup_file="/hosts/linuxhome/mutant9/tmp/sam/Prokaryotes/Caulobacter_crescentusR1/backups/backup01000000.txt";
// const string anctrace_file="/hosts/linuxhome/mutant9/tmp/sam/Prokaryotes/Caulobacter_crescentusR1/ancestors/anctrace01000000.txt";
//const string genome_file="";
//const string genestate_file="";
const string backup_file="";
const string anctrace_file="";

//Below is the stages as defined in the model version 2.1.
// const bool StageTargets[4][5] = {
//   false, false, true, false, true,      // 0 0 1 0 1    G1 -> S
//   false, true, false, false, false,     // 0 1 0 0 0    S -> G2
//   true, false, false, false, false,     // 1 0 0 0 0    G2 -> M
//   true, false, false, true, true        // 1 0 0 1 1    M -> G1
//                                         // 1-CtrA 2-GcrA 3-DnaA 4-CcrM 5-SciP
// };

//In version 2.2, they were redefined as below.
const bool StageTargets[4][5] = {
  true, false, false, true, true,       // 1 0 0 1 1    G1
  false, false, true, false, true,      // 0 0 1 0 1    S
  false, true, false, false, false,     // 0 1 0 0 0    G2
  true, false, false, false, false,     // 1 0 0 0 0    M
                                        // 1-CtrA 2-GcrA 3-DnaA 4-CcrM 5-SciP
};

#endif
