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
#include <boost/next_prior.hpp>
#include <cstring>
#include <sstream>
#include <fstream>
#include "dSFMT.h"	// Kirsten used it for some functions, so lets just stick with it
//#include "Xmgrace.hh"	// Grace pipeline added by Brem
#include <sys/time.h>
#include <ctime>
#include <typeinfo>
#include <string>

// #include "/home/sam/Programmes/nvwa-1.1/nvwa/debug_new.h"

#define toDigit(c) (c-'0')  // Converts char to digit

//Constants for Gene.hh
const int typeseq_length = 10;
const int binding_length = 20;

//Constants defined here, used in Population.hh
//Currently gives a Segmentation fault if NR and NC are not equal!
const int NR=50;  //While there is no diffusion, there seems to be no reason to make the field too big.
const int NC=550;

//Constants defined here, used in Genome.cc
const int init_nr_gene_types = 10;
const int init_nr_tfbs_per_gene = 3;
const int WeightRange = 3;  //TFBS weights can range from -WeightRange to +WeightRange.
const float repl_step_size = 80.;	//Effective rcs can be below 1, turning into a probability of replicating a single bead.
const bool replicate_entire_genes = false;
const bool repl_step_noise = false;
const int replication_time = 1;
const int nr_household_genes = 50;
const bool model_volume = false;

const bool relative_replication = false;
const int rel_repl_full = 65;	//If we're doing relative replication, how many nutrients are considered to be needed for replication of the entire genome.

const bool no_binding_noise = false;
const int binding_threshold = 17;	//If we're not doing binding noise, at what Hamming distance do genes start to bind. NOTE a possible alternative is to use the probabilities to scale all binding interactions, i.e. let everything always bind but proportional to its affinity (don't know if that is really necessary).

// Mutations
const bool type_mutations = false;
const double regulator_typeseq_mu = 0.00001;
const bool regtype[5][typeseq_length] =
{
	true, false, true, false, true, false, true, false, true, false,
	false, false, true, true, false, false, true, true, false, false,
	false, false, false, true, true, true, false, false, false, true,
	true, true, true, true, false, false, false, false, true, true,
	false, false, false, false, false, true, true, true, true, true
};

const double gene_threshold_mu = 0.0005;
const double gene_activity_mu = 0.0005;
const double gene_binding_domain_mu = 0.0001;

const double tfbs_binding_site_mu = 0.0001;
const double tfbs_activity_mu = 0.0005;

const double gene_duplication_mu = 0.0005;
const double gene_deletion_mu = 0.0005;
const double gene_innovation_mu = 0.0005; //I set this 10x lower than other mutation rates on purpose.
const double gene_destruction_mu = 0.0000;  //The opposite of an innovation; delete random gene from the genome. Note that this includes its binding sites. If you want to balance gene innovation mutations (which exclude binding sites), set the bead_destruction_mu to the same rate as gene, tfbs and house innovation combined.
const double gene_shuffle_mu = 0.001;

const double tfbs_duplication_mu = 0.0005;
const double tfbs_deletion_mu = 0.0005;
const double tfbs_innovation_mu = 0.005;
const double tfbs_shuffle_mu = 0.001;

const double house_duplication_mu = 0.0001;
const double house_deletion_mu = 0.0001;
const double house_innovation_mu = 0.000;
const double house_shuffle_mu = 0.001;

const double bead_destruction_mu = 0.0; //Opposite of all innovations combined (see above) -  on average there is no net increase or decrease of genome size. However, this would result in a bias within the genome fractions of genes, houses, and TFBSs (due to different innovation rates).

const double k_zero = 0.0000001;
const double epsilon = 1.00;

//constants for Population.cc
const int TimeZero=0;
const int default_SimTime=2000000;
const int TimeTerminalOutput = 100;  //Note that this is also the check for extinction, so let it check!
const int TimeSaveGrid = 100; //How many timesteps to save the whole grid as raw data.
const int TimeSaveBackup = 10000;
const int TimePruneFossils = 1000;
const int TimeOutputFossils = 10000;

const double diffusion_rate = 0.;  // >1: multiple Margolus steps per time step, <1: probability of single Margolus step each time step.
const double death_rate = 0.02;
const double repl_rate = 1.0;
const int replication_neighbourhood = 3;  //i.e. a 3x3 grid represents the neighbourhood.
const int generation_sample = 10000;

const bool environmental_noise = false;  //Actually used in Genome.cc and Population.cc (for initialisation).
const double environmental_change_rate = 0.01;
const int environmental_variation = 10;

const bool resource_dependent_replication = true;
const int env_blocks = 11;
const double gradient[11] = {0., 10., 20., 30., 40., 50., 60., 70., 72., 75., 78.};
const int stats_in_blocks = 11;		//In how many blocks should the field be split up to collect stats. Useful for analyses; for normal simulation it would be best to keep it at "1" (regard field as single block for stats).
// const int env_blocks = 1;
// const double gradient[1] = {50.};
// const int stats_in_blocks = 1;

//Protocol for division:
// 0 - overgrow neighbour.
// 1 - compete with neighbour.
// 2 - wait for empty site.
const int DivisionProtocol = 2;

//Penalties for "wrong" cell-cycles:
// 0 - no penalty but also not allowed to keep replicating or attempting mitosis if in "S" or "M".
// 1 - abortion; cell reset to Stage 1 and replicated beads are lost.
// 2 - death upon division; cell is marked for death when it attempts to divide.
// 3 - immediate death; cell is marked for immediate death.
const int ShortReplProtocol = 0;
const int EarlyMitProtocol = 3;
const int BadUpdProtocol = 0;

// Protocol settings: DivisionProtocol, ShortReplProtocol, EarlyMitProtocol, BadUpdProtocol.
// Muts_At_Divo:    2, 0, 0, 0
// Risk_Wait_Div:   2, 0, 2, 0
// Commit_Div:      0, 0, 2, 0
// Grow_Or_Stall:   1, ?, ?, 0

// C. pneumoniae:   ~0, 1, 2, 0   (I did not actually use 0 for DivisionProtocol, but "free waiting" in "M").
// C. trachomatis:  0, 0, 2, 0

//constants for Prokaryote.cc

//Variables defined in World.cc
extern int Time;
extern int initial_seed;
extern unsigned long long seed_draws;
extern string folder;
extern bool mutational_neighbourhood;
extern bool mutational_scanpath;
extern bool attractor_landscape;
extern bool follow_single_individual;
extern int NrMutants;
extern int NrInitialStates;
extern string genome_init;
extern string genestate_init;
extern string backup_reboot;
extern string anctrace_reboot;
extern bool mutations_on;
extern double init_env;
extern int SimTime;

extern dsfmt_t dsfmt;
inline double uniform()
{
  seed_draws ++;
  return dsfmt_genrand_close_open(&dsfmt);
}

const string genome_file="";
const string genestate_file="";
const string backup_file="";
const string anctrace_file="";

//The current definition of the stages.
const bool StageTargets[4][5] = {
  true, false, false, true, true,       // 1 0 0 1 1    G1
  false, false, true, false, true,      // 0 0 1 0 1    S
  false, true, false, false, false,     // 0 1 0 0 0    G2
  true, false, false, false, false,     // 1 0 0 0 0    M
                                        // 1-CtrA 2-GcrA 3-DnaA 4-CcrM 5-SciP
};

//Below is the stages as defined in the earliest version of the model.
// const bool StageTargets[4][5] = {
//   false, false, true, false, true,      // 0 0 1 0 1    G1 -> S
//   false, true, false, false, false,     // 0 1 0 0 0    S -> G2
//   true, false, false, false, false,     // 1 0 0 0 0    G2 -> M
//   true, false, false, true, true        // 1 0 0 1 1    M -> G1
//                                         // 1-CtrA 2-GcrA 3-DnaA 4-CcrM 5-SciP
// };

#endif
