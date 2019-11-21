#include "Prokaryote.hh"
#include "stdio.h"
#include <stdlib.h>
#include <png.h>
#include <zlib.h>

Prokaryote::Prokaryote() {
}

Prokaryote::~Prokaryote() {
	if(G!=NULL){
		delete (G);
		G=NULL;
	}
}

void Prokaryote::InitialiseProkaryote(){
	EmptyProkaryote();
	if(genome_init != "")	G->ReadInitialGenome();
	else	G->InitialiseRandomGenome();
}

void Prokaryote::ClonePPFromPP(Prokaryote* PP_template){
	EmptyProkaryote();
	//Copy genome from template. This includes the expression levels of each gene.
	G->CloneGenome(PP_template->G);
	Stage = PP_template->Stage;
}

void Prokaryote::Replicate(Prokaryote* PP_parent){
	EmptyProkaryote();	//Currently this is quite an overkill as first the standard random genome is initiated only to be fully replaced by the copy from the parent genome
	//Copy parental genome. Includes expression levels, so I will in the future change this function into something like DivideGenome(), as the genome is already replicated.
	G->CloneGenome(PP_parent->G);
	//Do mutations.
	if(mutations_on) G->MutateGenome();		//Beware that this might change the type of genes, so that the expression should also be updated after a succesfull mutation.
	//Importantly, the parent also returns to the D-stage (0).
	PP_parent->Stage = 0;
	PP_parent->ready_for_replication = false;
}

void Prokaryote::EmptyProkaryote()
{
	//Create the genome of the prokaryote
	G = NULL;
	G = new Genome();
	//Starts at first stage again
	Stage = 0;
}

void Prokaryote::UpdateCellCycle()	//Check whether changes in GeneStates make us go forward in the cell cycle.
{
	int expression;
	Genome::gene_iter git;

	for (int s=0; s<4; s++)
	{
		if (s==Stage)
		{
			int match_next_state = 0;
			for (int g=0; g<5; g++)	//g is the type of gene that we are looking for, ie G0-G4; but these may be shuffled (or some missing) from the actual GeneStates.
			{
				git = find(G->GeneTypes->begin(), G->GeneTypes->end(), g);
				if (git == G->GeneTypes->end())	expression = 0;
				else
				{
					int index = distance(G->GeneTypes->begin(), git);
					expression = G->GeneStates->at(index);
				}

				if(  (expression==0 && !StageTargets[s][g])  ||  (expression!=0 && StageTargets[s][g])  )
				{
					match_next_state++;
				}
				else	break;
			}
			if (match_next_state==5)
			{
				Stage++;	//You have reached the next stage.
				break;		//One stage improvement is more than enough for one timestep :)
			}
			break;
		}
	}
}

double Prokaryote::Fitness()
{
	double fitness = 1.0;
	for(int i=0; i<5; i++)
	{
		if (G->GeneStates->at(i) != TargetExpression[i])	fitness -= 1./5.;
	}
	return pow(fitness,2);
}
