#ifndef PopulationHeader
#define PopulationHeader

#include "Prokaryote.hh"
#include "Genome.hh"
#include "Gene.hh"
#include "TFBS.hh"
#include "Header.hh"
#include "FossilRecord.hh"
#include <cstdio>
#include <stdlib.h>

class Population
{
	public:
		Prokaryote* PPSpace[NR][NC];
		FossilRecord* Fossils;

		int p_nr_proks_;	//Probably need this...

		typedef std::list<int>::iterator iter;
		typedef std::list<Prokaryote*>::iterator iterpps;

		unsigned long long p_id_count_;	// Counter for all agents

		//Variables for output.
		int nr_birth_events;
		int nr_first_births;
		int cum_time_alive;
		double cum_fit_def;

		//These are for looking at the occurrence of evolution.
		Prokaryote* OldGeneration[generation_sample];

		int Environment;	//Add noise for a given period of time.

		Population();
		~Population();

		void InitialisePopulation();
		void ContinuePopulationFromBackup();
		void ReadBackupFile();
		void ReadAncestorFile();

		void ReproduceMasterGenome();
		void FollowSingleIndividual();
		void ExploreAttractorLandscape();
		void UpdatePopulation();
		void DeathOfProkaryote(int i, int j);
		void SetEnvironment();

		void PruneFossilRecord();
		double MatrixDistance(Prokaryote* PP1, Prokaryote* PP2);
		void PrintFieldToFile();
		void PrintSampleToFile();
		void OutputBackup();
		void ShowGeneralProgress();
};

#endif
