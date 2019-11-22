#ifndef ProkaryoteHeader
#define ProkaryoteHeader

#include "Header.hh"
#include "Genome.hh"
#include "stdlib.h"

class Prokaryote{
	public:
		Genome* G;
		int Stage;
		bool ready_for_division;
		double fitness_deficit;

		typedef std::list<Prokaryote*>::iterator iteragent;

		Prokaryote();
		~Prokaryote();

		void InitialiseProkaryote();
		void ClonePPFromPP(Prokaryote* PP_template, int tot_prok_count);
		void Replicate();
		void Mitosis(Prokaryote* parent, int tot_prok_count);
		void EmptyProkaryote();
		void PrintData(bool include_genome_data);
		void UpdateCellCycle();
};

#endif
