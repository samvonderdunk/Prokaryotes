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

void Prokaryote::Replicate(int env, int res)
{
	if (G->pos_fork != G->pos_anti_ori)	//If the fork is has reached the opposite of ORI of the genome, there is nothing to replicate.
	{
		if (uniform() <= 1.0)	//later the chance that replication proceeds one step depends on several things.
		{
			G->ReplicateGenomeStep(env, res);
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

	if (house_duplication_mu > 0.0 || house_deletion_mu > 0.0)	//We only have to check the number of household genes if they can actually
	{
		fitness_deficit = abs(nr_household_genes - G->gnr_houses) / (float)10;
	}
}

void Prokaryote::Abortion()
{
	G->AbortChildGenome();

	Stage = 0;
	ready_for_division = false;
	time_replicated = 0;
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
	time_replicated = 0;
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
	int s, match_next_state;

	for (s=0; s<5; s++)
	{
		if (s==Stage)
		{
			// if (Stage == 4 || (Stage == 2 && G->pos_fork != G->pos_anti_ori))	//If you were in Stage 4, you have to earn it again; you are set back to stage 3 and evaluated for matching stage 4 anew.
			if (Stage == 4 || (Stage == 2 && (time_replicated < replication_time || G->pos_fork != G->pos_anti_ori)))	//Substitute for the above line to learn prokaryotes the trick of S-stage extension.
			// if(Stage == 2 && (time_replicated < replication_time || G->pos_fork != G->pos_anti_ori))	//Prokaryotesv2.3: You don't stay in M-stage anymore, so the above lines are deprecated.
			{
				Stage--;
				s--;	//This means we will evaluate a Stage-4 cell for Stage 4 again.
			}

			match_next_state = G->MatchNextState(s);

			if (match_next_state==5)
			{
				Stage++;	//You have reached the next stage.
			}

			// else	//If you don't reach your normal next-state, we check whether you have reached the M-stage. Then you are updated to M-stage immediately. If you do not go there from Stage 3 (G2, as you should) we put your time_replicated to 0. Thus, reaching M-stage too soon will always kill you. You are not allowed to skip any of G1, S (multiple steps) and G2.
			// {
			// 	match_next_state = G->MatchNextState(3);	//Evaluate for M-stage.
			// 	if (match_next_state==5)
			// 	{
			// 		time_replicated = 0;	//The death penalty.
			// 		Stage = 4;
			// 	}
			// }

			/*	//Use this code to penalize not matching your next state and if you want, to set expression to the correct state.
			else	//Let's try to penalize everything that does not fit into our cell-cycle scheme: G1-S-S-S-S-(etc.)-G2-M-M-M-M-(mitosis)-G1-S-S-S-S...
			{
				fitness_deficit += 0.01 * (5-match_next_state);
				Stage = s+1;

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
						gene->expression = (StageTargets[s][gene->type-1]==true) ? 1:0;	//Set this gene to what it should be in S-stage.
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
			*/

			break;		//One stage evaluation is more than enough for one timestep :)
		}
	}
}
