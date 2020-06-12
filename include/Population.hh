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
		typedef std::list<unsigned long long>::iterator iterull;
		typedef std::list<Prokaryote*>::iterator iterpps;
		typedef std::pair<int,int> coords;	//Allow to define row, col pair of int's (use for functions).

		unsigned long long p_id_count_;	// Counter for all agents

		//Variables for output.
		int nr_birth_events[stats_in_blocks];
		int nr_first_births[stats_in_blocks];
		int cum_time_alive[stats_in_blocks];
		double cum_fit_def[stats_in_blocks];
		int nr_death_cycles[stats_in_blocks];

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
		void MargolusDiffusion();
		bool IsReadyToDivide(int i, int j, int nrow, int ncol);
		void DeathOfProkaryote(int i, int j);

		void NoiseEnvironment();
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
		void PrintStatPerBlock(string stat_name, int* stat_array, int* scaler_array);
		void PrintStatPerBlock(string stat_name, const double* stat_array, int* scaler_array);	//Version that accepts a double array.
};

#endif
