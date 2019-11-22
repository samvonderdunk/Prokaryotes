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

void Prokaryote::PrintData(bool include_genome_data)
{
	printf("##############################################\n");
	printf("Prokaryote #%llu:\n", fossil_id);
	if(Ancestor==NULL)	printf("Generation 0\n");
	else	printf("Child of #%llu\n", Ancestor->fossil_id);
	printf("Stage = %d\n", Stage);
	printf("Time of birth = %d\n", time_of_appearance);
	printf("Fitness deficit = %f\n", fitness_deficit);
	printf("It is %s.\n", (alive)? "alive":"dead");
	printf("It is %sready for division.\n", (ready_for_division)? "":"not ");
	printf("It is %sa mutant.\n", ((mutant)? "":"not " ));
	printf("It is %sbearing a mutant child.\n", ((mutant_child)? "":"not " ));
	printf("It will %sbe saved in the graveyard.\n", ((saved_in_graveyard)? "":"not "));
	printf("----------------------------------------------\n");
	if(include_genome_data)
	{
		printf("Genome:\n%s\n", G->PrintContent(NULL, false, false).c_str());
		printf("GeneStates:\n%s\n", G->PrintGeneStateContent().c_str());
		printf("GeneTypes:\n%s\n", G->PrintGeneTypeContent().c_str());
	}
	printf("##############################################\n");
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
			for (int g=1; g<=5; g++)	//g is the type of gene that we are looking for, ie G1-G5; but these may be shuffled (or some missing) from the actual GeneStates.
			{
				git = find(G->GeneTypes->begin(), G->GeneTypes->end(), g);
				if (git == G->GeneTypes->end())	expression = 0;
				else
				{
					int index = distance(G->GeneTypes->begin(), git);
					expression = G->GeneStates->at(index);
				}

				if(  (expression==0 && !StageTargets[s][g-1])  ||  (expression!=0 && StageTargets[s][g-1])  )
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
