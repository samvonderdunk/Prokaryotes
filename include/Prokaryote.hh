#ifndef ProkaryoteHeader
#define ProkaryoteHeader

#include "Header.hh"
#include "Genome.hh"
#include "stdlib.h"

class Prokaryote{
	public:
		Genome* G;
		int Stage;
		bool ready_for_replication;

		typedef std::list<Prokaryote*>::iterator iteragent;

		Prokaryote();
		~Prokaryote();

		void InitialiseProkaryote();
		void ClonePPFromPP(Prokaryote* PP_template);
		void Replicate(Prokaryote* PP_parent);
		void EmptyProkaryote();

		void UpdateCellCycle();
		double Fitness();
};

#endif
