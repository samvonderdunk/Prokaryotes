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
	mutant = true;	//The first will be deemed a mutant so that the ancestor trace always comes back to at least one individual of the initial batch.
	if(genome_init != "")	G->ReadInitialGenome();
	else	G->InitialiseRandomGenome();
}

void Prokaryote::ClonePPFromPP(Prokaryote* PP_template, int tot_prok_count)
{
	EmptyProkaryote();
	//Copy genome from template. This includes the expression levels of each gene.
	G->CloneGenome(PP_template->G);
	Stage = PP_template->Stage;
	fossil_id = tot_prok_count;
	time_of_appearance = Time;
	mutant = PP_template->mutant;	//If you clone a mutant (or the first prokaryote), its clone will also count as a mutant.
}

void Prokaryote::Replicate(int env)
{
	if (G->pos_fork != G->pos_anti_ori)	//If the fork is has reached the opposite of ORI of the genome, there is nothing to replicate.
	{
		if (uniform() <= 1.0)	//later the chance that replication proceeds one step depends on several things.
		{
			G->ReplicateGenomeStep(env);
		}
	}
}

void Prokaryote::Mitosis(Prokaryote* parent, int tot_prok_count)
{
	EmptyProkaryote();
	fossil_id = tot_prok_count;
	time_of_appearance = Time;

	G->SplitGenome(parent->G);

	if (G->mutant_genome)	mutant = true;
	else	mutant = false;

	if (parent->mutant)	Ancestor = parent;	//If your parent was a mutant (i.e. its genome holds a mutation with respect to its parent), then your immediate ancestor is your parent.
	else	Ancestor = parent->Ancestor;	//If your parent was not a mutant (its genome is the same as its parent), then point to its MRCA.

	parent->Stage = 0;
	parent->ready_for_division = false;
	// parent->fitness_deficit = 0.;	//It can try to replicate better next time.
	parent->nr_offspring++;
	parent->time_replicated = 0;
}

void Prokaryote::EmptyProkaryote()
{
	//Create the genome of the prokaryote
	G = NULL;
	G = new Genome();
	//Starts at first stage again
	Stage = 0;
	ready_for_division = false;
	nr_offspring = 0;
	fossil_id = 0;
	time_of_appearance = 0;
	Ancestor = NULL;
	mutant = false;
	alive = true;
	fitness_deficit = 0.;
	saved_in_graveyard = false;
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
	printf("It will %sbe saved in the graveyard.\n", ((saved_in_graveyard)? "":"not "));
	printf("----------------------------------------------\n");
	if(include_genome_data)
	{
		printf("Genome:\n%s\n", G->PrintContent(NULL, true, false).c_str());
		printf("GeneStates:\n%s\n", G->PrintGeneStateContent(false).c_str());
		printf("GeneTypes:\n%s\n", G->PrintGeneTypeContent().c_str());
	}
	printf("##############################################\n");
}

void Prokaryote::UpdateCellCycle()	//Check whether changes in GeneStates make us go forward in the cell cycle.
{
	int expression, index;
	Genome::gene_iter git;
	int st;

	for (int s=0; s<5; s++)
	{
		if (s==Stage)
		{
			if (Stage == 4)	//If you were in Stage 4, you have to earn it again; you are set back to stage 3 and evaluated for matching stage 4 anew.
			{
				Stage = 3;
				s--;	//This means we will evaluate a Stage-4 cell for Stage 4 again.
			}
			int match_next_state = 0;
			for (int g=1; g<=5; g++)	//g is the type of gene that we are looking for, ie G1-G5; but these may be shuffled (or some missing) from the actual GeneStates.
			{
				git = find(G->GeneTypes->begin(), G->GeneTypes->end(), g);
				if (git == G->GeneTypes->end())	expression = 0;
				else
				{
					index = distance(G->GeneTypes->begin(), git);
					expression = G->GeneStates->at(index);
				}

				if(  (expression==0 && !StageTargets[s][g-1])  ||  (expression!=0 && StageTargets[s][g-1])  )
				{
					match_next_state++;
				}
				else	break;
			}
			if (match_next_state==5)	//In principle you develop if you match your next state (and no worries if you don't immediately). When you are in M-stage, however, you should match M-stage, or you get forced back into M-stage with a fitness cost.
			{
				if ((Stage == 2 && G->pos_fork != G->pos_anti_ori))	//Check that you are done replicating.
				{
					fitness_deficit += 0.1;

					Stage--;	//Keep cell in S- or M-stage virtually.

					git = G->GeneTypes->begin();	//Remove expression of cell-cycle genes.
					while (git != G->GeneTypes->end())
					{
						if (*git < 6)
						{
							index = distance(G->GeneTypes->begin(), git);
							G->GeneStates->at(index) = 0;
						}
						git++;
					}

					Genome::iter it = G->BeadList->begin();	//Now keep it in there physically.
					while (it != G->BeadList->end())
					{
						if(G->IsGene(*it) && (*it)->type < 6)
						{
							Gene* gene = dynamic_cast<Gene*>(*it);
							gene->expression = (StageTargets[s-1][gene->type-1]==true) ? 1:0;	//Set this gene to what it should be in S-stage.
							git = find(G->GeneTypes->begin(), G->GeneTypes->end(), gene->type);
							if(git != G->GeneTypes->end())
							{
								index = distance(G->GeneTypes->begin(), git);
								G->GeneStates->at(index) += gene->expression;
							}
						}
						it++;
					}
				}

				Stage++;	//You have reached the next stage.

				if (Stage == 3)	//If you opt out of the S-stage you better make sure that everything had a chance to mutate.
				{
					Genome::mut_iter mit = G->MutationList->begin();
					while (mit != G->MutationList->end())
					{
						assert((*mit));
						mit++;
					}

					delete G->MutationList;
					G->MutationList = NULL;
				}
				break;		//One stage improvement is more than enough for one timestep :)
			}
			break;
		}
	}
}
