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

		int p_id_count_;	// Counter for all agents

		Population();
		~Population();

		void InitialisePopulation();
		void ReproduceMasterGenome();
		void UpdatePopulation();
		void PruneFossilRecord();
		void PrintFieldToFile();
		void PrintSampleToFile();
		void ShowGeneralProgress();
};

#endif
