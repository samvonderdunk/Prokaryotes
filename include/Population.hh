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
		typedef std::pair<int,int> coords;	//Allow to define row, col pair of int's (use for functions).

		unsigned long long p_id_count_;	// Counter for all agents

		//Variables for output.
		int nr_birth_events;
		int nr_first_births;
		int cum_time_alive;
		double cum_fit_def;

		//These are for looking at the occurrence of evolution.
		Prokaryote* OldGeneration[generation_sample];

		double Environment;	//Add noise for a given period of time.

		Population();
		~Population();

		void InitialisePopulation();
		void ContinuePopulationFromBackup();
		void ReadBackupFile();
		void ReadAncestorFile();

		void ReproduceMasterGenome();
		void ScanMutationalPath();
		bool CompareExpressionProgression(Prokaryote* PP1, Prokaryote* PP2, int time_window);
		void FollowSingleIndividual();
		void ExploreAttractorLandscape();

		void UpdatePopulation();
		bool IsReadyToDivide(int i, int j, int nrow, int ncol);
		void DeathOfProkaryote(int i, int j);

		void SetEnvironment();
		void GradientEnvironment(int i, int j);
		double CollectResource(int i, int j, double Environment);
		coords PickNeighbour(int i, int j);
		int NeighbourhoodDensity(int i, int j);
		void ResetProgressCounters();

		void PruneFossilRecord();

		double MatrixDistance(Prokaryote* PP1, Prokaryote* PP2);
		bool QualitativeNetworkChange(Prokaryote* PP1, Prokaryote* PP2);
		void PrintFieldToFile();
		void PrintSampleToFile();	//Use this for samples of the field.
		void OutputBackup();
		void ShowGeneralProgress();
};

#endif
