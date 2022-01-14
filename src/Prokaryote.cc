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
	time_stationary = 0;
	priviliges = false;
	maturing_time = 0;
}

void Prokaryote::InitialiseProkaryote(){
	EmptyProkaryote();
	mutant = true;	//The first will be deemed a mutant so that the ancestor trace always comes back to at least one individual of the initial batch.
	if(genome_init != "")	G->ReadInitialGenome();
	else	G->InitialiseRandomGenome();
}

void Prokaryote::ClonePPFromPP(Prokaryote* PP_template, unsigned long long tot_prok_count)
{
	EmptyProkaryote();
	//Copy genome from template. This includes the expression levels of each gene.
	G->CloneGenome(PP_template->G);
	Stage = PP_template->Stage;
	fossil_id = tot_prok_count;
	time_of_appearance = Time;
	mutant = PP_template->mutant;	//If you clone a mutant (or the first prokaryote), its clone will also count as a mutant.
}

void Prokaryote::Replicate(double resource)
{
	if (G->pos_fork != G->pos_anti_ori)	//If the fork is has reached the opposite of ORI of the genome, there is nothing to replicate.
	{
		G->ReplicateGenomeStep(resource);
		time_replicated++;
	}
}

void Prokaryote::Mitosis(Prokaryote* parent, unsigned long long tot_prok_count)
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
		// fitness_deficit = min(1., (nr_household_genes - G->gnr_houses) / (double)10);
	}
}

void Prokaryote::UpdateCellCycle()	//Check whether changes in GeneStates make us go forward in the cell cycle.
{
	int evaluate_stage = Stage;	//Note that Stage includes "D", so "S" corresponds to Stage 2 but to EvaluateStage 1.
	priviliges = true;	//Priviliges are only removed for one turn that you did not acquire the desired expression...

		//M and S expression has to be maintained actively (in protocols where division does not require an empty site, individuals will never get here in "M" expression; but otherwise staying in "M" is not free).
		//Stage 4 and 2 cells will be evaluated for those stages again.
	if (Stage == 4 || (Stage == 2 && (time_replicated < replication_time || G->pos_fork != G->pos_anti_ori)))		evaluate_stage--;

	/*Evaluation*/

            //You have reached the next stage.
	if (G->MatchNextState(evaluate_stage) == 5)    Stage = evaluate_stage + 1;
	else
	{
			//The worst penalty should come first...
		if (G->MatchNextState(3) == 5)    UpdatePenalty(EarlyMitProtocol);
			//The strictest evaluation criterium first (i.e. NOT staying in "S" when you're not done replicating).
		else if (Stage == 2 && (time_replicated < replication_time || G->pos_fork != G->pos_anti_ori))    UpdatePenalty(ShortReplProtocol);

		else  UpdatePenalty(BadUpdProtocol);

			//Even if we don't do competition upon division it does not hurt to track time_stationary by default.
		if (G->MatchNextState(0) == 5 || G->MatchNextState(2) == 5)   time_stationary++;
	}
}

void Prokaryote::UpdatePenalty(int protocol)
{
	//Penalties sorted from the leanest to the strictest.
	if (protocol == 0)				priviliges = false;	//No penalty, but does not allow you to keep replicating or attempting in Stage 2 or 4. Since the next penalties are all worse, we don't have to set the priviliges to false in those cases.
	else if (protocol == 1)		Abortion();
	else if (protocol == 2)		Stage = 5;	//Sentenced to wait until mitosis and then die (it won't update anymore). NOTE that priviliges are not turned off, because that would prevent the cell from attempting division!
	else if (protocol == 3)		Stage = 6;	//Sentenced to death immediately (see Population.cc).
}

void Prokaryote::Abortion()
{
	G->AbortChildGenome();

	Stage = 0;
	ready_for_division = false;
	time_replicated = 0;
	time_stationary = 0;
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
