#include "Header.hh"
#include "Population.hh"
#include "Prokaryote.hh"
#include "Genome.hh"
#include "Gene.hh"
#include "TFBS.hh"
#include "dSFMT.h"

dsfmt_t dsfmt;
int Time;
int initial_seed = time(0);
string folder = "/linuxhome/tmp/sam/Prokaryotes/";
bool mutational_neighbourhood = false;
bool mutational_scanpath = false;
bool attractor_landscape = false;
bool follow_single_individual = false;
int NrMutants = 0;
int NrInitialStates = 0;
string genome_init = genome_file;
string genestate_init = genestate_file;
string backup_reboot = backup_file;
string anctrace_reboot = anctrace_file;
bool mutations_on = true;
int init_env = 0;
int SimTime = default_SimTime;

void Setup(int argc, char** argv);

int main(int argc, char** argv) {

	/* ############## Setup ############## */
	printf("\n\033[93m### Setup ###\033[0m\n");
	Population* P;
	Setup(argc, argv);
	printf("Function call: ");
	for(int q=0; q<argc; q++)	printf("%s ", argv[q]);
	dsfmt_init_gen_rand(&dsfmt, initial_seed);	//Used to seed uniform().
	srand(initial_seed);	//Used to seed random_shuffle(...).
	printf("\b\nSetup completed...\n\n");

	/* ############## Simulate Mutants ############## */
	if (mutational_neighbourhood)
	{
		printf("\033[93m### Start ###\033[0m\n");
		P = new Population();
		P->ReproduceMasterGenome();
	}

	else if (mutational_scanpath)
	{
		printf("\033[93m### Start ###\033[0m\n");
		P = new Population();
		P->ScanMutationalPath();
	}

	else if (attractor_landscape)
	{
		printf("\033[93m### Start ###\033[0m\n");
		P = new Population();
		P->ExploreAttractorLandscape();
	}

	else if (follow_single_individual)
	{
		printf("\033[93m### Start ###\033[0m\n");
		P = new Population();
		P->FollowSingleIndividual();
	}

	else
	{
		/* ############## Initialisation ############## */
		printf("\033[93m### Initialisation ###\033[0m\n");
		P = new Population();
		if(backup_reboot != "")	P->ContinuePopulationFromBackup();
		else	P->InitialisePopulation();
		printf("Initialisation completed...\n\n");

		/* ############## Simulation ############## */

		printf("\033[93m### Simulation ###\033[0m\n");
		for(Time=TimeZero; Time<SimTime+1; Time++){	//We do one extra step, because output is generated at the beginning of a step, such that time=0 is the field as it is initialised.
			P->UpdatePopulation();		//Main next-state function, updating the population.
		}
		// cout << P->nr_birth_events << endl;
		//Make sure that you save all possible things in the last timestep, if you did not already choose your parameters such.
		Time--;
		if(SimTime%TimeTerminalOutput!=0)	P->ShowGeneralProgress();
		if(SimTime%TimeSaveGrid!=0)
		{
			if(NR*NC > 3000)	P->PrintSampleToFile();
			else	P->PrintFieldToFile();
		}
		if(SimTime%TimeSaveBackup!=0)	P->OutputBackup();
		if(SimTime%TimePruneFossils!=0)	P->PruneFossilRecord();
		if(SimTime%TimeOutputFossils!=0)	P->Fossils->ExhibitFossils();
		printf("Simulation completed...\n\n");
	}

	/* ############## End ############## */

	printf("\033[93m### End ###\033[0m\n");
	delete P;
	P = NULL;
	printf("Prokaryotes completed...\n\n");

}



void Setup(int argc, char** argv) {

	string ReadOut, command;
	bool project_name_found = false;
	bool initial_seed_set = false;

	for(int i=1;i<argc;i++)	//Loop through input arguments.
	{
		ReadOut = (char*) argv[i];	//There does not seem to be a quicker way to compare the input arguments with a string.

		//Let user define initial_seed. Otherwise defaults to time(0).
		if(ReadOut=="-s" && (i+1)!=argc)
		{
			initial_seed = atoi(argv[i+1]);
			initial_seed_set = true;
			printf("Seed = %i\n", initial_seed);
			i++;
			continue;
		}

		//Let user define subdirectory for the project
		else if(ReadOut=="-p" && (i+1)!=argc)
		{
			folder += argv[i+1];
			project_name_found = true;
			i++;
			continue;
		}

		//If user wants to look at mutational neighbourhood, no simulation will be started, and a bunch of offspring is generated. Make sure to provide a genome file, otherwise the programme will use a randomly generated genome to seed the offspring.
		else if(ReadOut=="-M" && (i+1)!=argc)
		{
			mutational_neighbourhood = true;
			NrMutants = atoi(argv[i+1]);
			printf("Exploring mutational neighbourhood by generating %d children.\n", NrMutants);
			i++;
			continue;
		}

		else if(ReadOut=="-MS")
		{
			mutational_scanpath = true;
			printf("Walking along the neutral mutational path.\n");
		}

		//If user wants to look at attractor landscape, all possible states will be simulated for one timestep (or a sample of all possible states) and their resulting states will be given as output.
		else if(ReadOut=="-A" && (i+1)!=argc)
		{
			attractor_landscape = true;
			NrInitialStates = atoi(argv[i+1]);
			printf("Exploring attractor landscape by simulating %d initial states.\n", NrInitialStates);
			i++;
			continue;
		}

		//Follow a single immortal individual for many time steps.
		else if(ReadOut=="-S")
		{
			follow_single_individual = true;
			printf("Following a single, immortal individual through time.\n");
		}

		//Don't do mutations, i.e. no evolution.
		else if(ReadOut=="-nomut")
		{
			mutations_on = false;
			printf("Simulating without mutations.\n");
		}

		//Let user define input genome file on the command line (allows you to work together with snakemake at kindergarten).
		else if(ReadOut=="-i" && (i+1)!=argc)
		{
			genome_init = argv[i+1];
			printf("Genome input: %s\n", genome_init.c_str());
			i++;
			continue;
		}

		//Let user define input genestate file on the command line (again handy for snakemake). You can either define the path to a genestate file or input a letter corresponding to one of the four cell-cycle stages (G1, S, G2 or M).
		else if(ReadOut=="-g" && (i+1)!=argc)
		{
			genestate_init = argv[i+1];
			printf("Genestate input: %s\n", genestate_init.c_str());
			i++;
			continue;
		}

		else if(ReadOut=="-b" && (i+1)!=argc)
		{
			backup_reboot = argv[i+1];
			printf("Backup-file input: %s\n", backup_reboot.c_str());
			i++;
			continue;
		}

		else if(ReadOut=="-a" && (i+1)!=argc)
		{
			anctrace_reboot = argv[i+1];
			printf("Anctrace input: %s\n", anctrace_reboot.c_str());
			i++;
			continue;
		}

		else if(ReadOut=="-e" && (i+1)!=argc)
		{
			init_env = atoi(argv[i+1]);
			printf("Environment input: %d\n", init_env);
			i++;
			continue;
		}

		else if(ReadOut=="-t" && (i+1)!=argc)
		{
			SimTime = atoi(argv[i+1]);
			printf("Simulation time: %d\n", SimTime);
			i++;
			continue;
		}

		else	//Print usage/help.
		{
			printf("\n\033[93m### Prokaryotes --- usage ###\033[0m\nArgument options:\n   -p [project title]\t\tDefines folder for local storage\n   -s [seed]\t\t\tSet seed for random number generator (e.g. 211)\n   -i [initial genome]\t\te.g. MRCA.g\n   -g [initial expression]\te.g. MRCA_GS.g\n   -e [env]\t\t\tInitial environment (e.g. -3)\n   -b [backup file]\t\tStart from backup (e.g. /path/backup00090000.txt)\n   -a [ancestor file]\t\tContinue ancestor trace (e.g. /path/anctrace00090000.txt)\n   -nomut\t\t\tNo mutations\n   -t [max. time]\t\tSet simulation time (e.g. 100)\n Programmes:\n   -M [nr_mutants]\tGenerate mutants\n   -MS\t\t\tScan neutral mutational path\n   -A [nr_states]\tSimulate state-space transitions [nr. of initial states = max(nr_states, total nr. unique states)]\n   -S\t\t\tFollow single immortal individual/lineage through time [simulating until Time==SimTime]\n");
			exit(1);
		}
	}

	if (!initial_seed_set)	printf("Seed = %li\n", time(0));
	if (!project_name_found)	folder += "Project_Name";	//I did not manage to give the date as an extension to the folder.

	command = "mkdir -p " + folder;
	system(command.c_str());
	printf("Folder = %s\n", folder.c_str());
	//Automatically set up a subdirectory for snapshots of the grid (not images but raw data).
	command = "mkdir -p " + folder + "/snapsamples";
	system(command.c_str());
	command = "mkdir -p " + folder + "/backups";
	system(command.c_str());
	command = "mkdir -p " + folder + "/ancestors";
	system(command.c_str());

}
