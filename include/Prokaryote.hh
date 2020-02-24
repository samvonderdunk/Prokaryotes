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

		int time_replicated;

		//For output.
		int nr_offspring;	//We can measure the length of our beasts' life cycles (i.e. Time - time_of_appearance for all individuals at the point of doing Mitosis with the same nr_offspring).

		//For the fossil record.
		int time_of_appearance;
		unsigned long long fossil_id;	//Now it should be 64-bit (32-bit unsigned would already be about 200 times as big as the last id you get out of a 100x100 simulation of 100k AUT).
		Prokaryote* Ancestor;
		bool mutant;
		bool alive;
		bool saved_in_graveyard;	//For those individuals who are kept around a bit to calculate MSD.

		typedef std::list<Prokaryote*>::iterator iterpps;

		Prokaryote();
		~Prokaryote();

		void InitialiseProkaryote();
		void ClonePPFromPP(Prokaryote* PP_template, int tot_prok_count);
		void Replicate(int env, int res);
		void Mitosis(Prokaryote* parent, int tot_prok_count);
		void Abortion();
		void EmptyProkaryote();
		void PrintData(bool include_genome_data);
		void UpdateCellCycle();
};

#endif
