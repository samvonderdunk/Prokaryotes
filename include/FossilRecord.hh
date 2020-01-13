#ifndef FossilHeader
#define FossilHeader

#include "Prokaryote.hh"
#include "Genome.hh"
#include "Gene.hh"
#include "TFBS.hh"
#include "Header.hh"
#include <cstdio>
#include <stdlib.h>

class FossilRecord
{
	public:
		std::list<Prokaryote*> FossilList;

		typedef std::list<Prokaryote*>::iterator iterpps;
  	FossilRecord();
  	~FossilRecord();
		void EraseFossil(unsigned long long fossilID);
  	void BuryFossil(Prokaryote *P);
		void ExhibitFossils();
};
#endif
