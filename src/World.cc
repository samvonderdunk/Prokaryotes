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
int TargetExpression[5] = {};
string folder = "/linuxhome/tmp/sam/Prokaryotes/";
bool mutational_neighbourhood = false;
int NrMutants = 0;
string genome_init = genome_file;
string genestate_init = genestate_file;
bool mutations_on = true;

void Setup(int argc, char** argv);

int main(int argc, char** argv) {

	/* ############## Setup ############## */

	printf("\n\033[93m### Setup ###\033[0m\n");
	Population* P;
	Setup(argc, argv);
	dsfmt_init_gen_rand(&dsfmt, initial_seed);	//Used to seed uniform().
	srand(initial_seed);	//Used to seed random_shuffle(...).
	printf("Setup completed...\n\n");

	/* ############## Simulate Mutants ############## */
	if (mutational_neighbourhood == true)
	{
		printf("\033[93m### Start ###\033[0m\n");
		P = new Population();
		P->ReproduceMasterGenome();
	}

	else
	{
		/* ############## Initialisation ############## */
		printf("\033[93m### Initialisation ###\033[0m\n");
		P = new Population();
		P->InitialisePopulation();
		printf("Initialisation completed...\n\n");

		/* ############## Simulation ############## */

		printf("\033[93m### Simulation ###\033[0m\n");
		for(Time=0; Time<SimTime+1; Time++){	//We do one extra step, because output is generated at the beginning of a step, such that time=0 is the field as it is initialised.
			P->UpdatePopulation();		//Main next-state function, updating the population.

		}

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

	// string genome_init = genome_file;
	// printf("Genome init: %s", genome_init.c_str());
	// string genestate_init = genestate_file;
	// printf("Genestate init: %s", genestate_init.c_str());


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
		if(ReadOut=="-p" && (i+1)!=argc)
		{
			folder += argv[i+1];
			project_name_found = true;
			i++;
			continue;
		}

		//If user wants to look at mutational neighbourhood, no simulation will be started, and a bunch of offspring is generated. Make sure to provide a genome file, otherwise the programme will use a randomly generated genome to seed the offspring.
		if(ReadOut=="-M" && (i+1)!=argc)
		{
			mutational_neighbourhood = true;
			NrMutants = atoi(argv[i+1]);
			printf("Exploring mutational neighbourhood by generating %d children.\n", NrMutants);
			i++;
			continue;
		}

		//Don't do mutations, i.e. no evolution.
		if(ReadOut=="-nomut")
		{
			mutations_on = false;
			printf("Simulating without mutations.\n");
		}

		//Let user define input genome file on the command line (allows you to work together with snakemake at kindergarten).
		if(ReadOut=="-i" && (i+1)!=argc)
		{
			genome_init = argv[i+1];
			printf("Genome input: %s\n", genome_init.c_str());
		}

		//Let user define input genestate file on the command line (again handy for snakemake). You can either define the path to a genestate file or input a letter corresponding to one of the four cell-cycle stages (G1, S, G2 or M).
		if(ReadOut=="-g" && (i+1)!=argc)
		{
			genestate_init = argv[i+1];
			printf("Genestate input: %s\n", genestate_init.c_str());
		}

	}

	if (!initial_seed_set)	printf("Seed = %li\n", time(0));
	if (!project_name_found)	folder += "Project_Name";	//I did not manage to give the date as an extension to the folder.

	command = "mkdir -p " + folder;
	system(command.c_str());
	printf("Folder = %s\n", folder.c_str());
	//Automatically set up a subdirectory for snapshots of the grid (not images but raw data).
	command = "mkdir -p " + folder + "/snapgrids";
	system(command.c_str());
	command = "mkdir -p " + folder + "/snapsamples";
	system(command.c_str());

}
