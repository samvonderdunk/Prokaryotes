#include "Genome.hh"

Genome::Genome() {
	BeadList=NULL;
	GeneStates=NULL;
	GeneTypes=NULL;
	gnr_genes=0;
	gnr_houses=0;
	g_length=0;
	pos_fork=0;
	pos_anti_ori=0;
	mutant_genome=false;
}

Genome::~Genome() {
	iter i;

	if(BeadList!=NULL) {
		i=BeadList->begin();
		while(i!=BeadList->end()) {
			delete (*i);
			i++;
		}

		i=BeadList->erase(BeadList->begin(),BeadList->end());
		delete BeadList;
		BeadList=NULL;
	}

	delete GeneStates;
	delete GeneTypes;
	GeneStates=NULL;
	GeneTypes=NULL;
}

void Genome::CopyPartOfGenome(iter begin, iter end)
{
	iter ii;
	Bead* bead;
	ii=begin;
	while(ii!=end)
	{
		bead=(*ii)->Clone();
		(*BeadList).push_back(bead);
		g_length++;
		ii++;
	}
}

void Genome::CopyPartOfGenomeToTemplate(iter begin, iter end, list<Bead*>* template_beadlist)
{
	iter ii;
	Bead* bead;

	ii=begin;
	while(ii!=end)
	{
		bead=(*ii)->Clone();
		(*template_beadlist).push_back(bead);
		ii++;
	}
}


void Genome::CloneGenome(const Genome* G_template)	//Used to copy genome of prokaryote exactly
{
	BeadList=new list<Bead*>();
	CopyPartOfGenome(G_template->BeadList->begin(),G_template->BeadList->end());
	GeneStates=new vector<int>(G_template->GeneStates->size(),0);
	GeneTypes=new vector<int>(G_template->GeneStates->size(),0);
	for(int i=0; (size_t)i<G_template->GeneStates->size(); i++)
	{
		GeneStates->at(i) = G_template->GeneStates->at(i);
		GeneTypes->at(i) = G_template->GeneTypes->at(i);
	}

	gnr_genes = G_template->gnr_genes;
	gnr_houses = G_template->gnr_houses;
	g_length = G_template->g_length;
	pos_fork = G_template->pos_fork;
	pos_anti_ori = G_template->g_length;
}

void Genome::DevelopChildrenGenomes(Genome* G_replicated)	//Function gets iterators of parental genome, but copies it to child and then acts on variables of child genome.
{
	g_length = BeadList->size();
	GeneTypes = new vector<int>();
	GeneStates = new vector<int>();
	Gene* gene;
	int g_length_before_mut = g_length;
	vector<bool>* MutationList;
	int del_length, dup_length, index, type_abundance;
	int* pdup_length, * pdel_length;
	iter it;

	//Clean up the variables of the parental genome left behind.
	it = BeadList->begin();
	while (it != BeadList->end())
	{
		if(IsGene(*it))
		{
			type_abundance = G_replicated->CountTypeAbundance(abs((*it)->type));
			assert(type_abundance != 0);	//We have not done any mutations yet, so each gene should still have a complementary copy in the parent.
			G_replicated->DecrementExpressionOfType(it);
			G_replicated->gnr_genes--;
		}
		else if(IsHouse(*it))
		{
			G_replicated->gnr_houses--;
		}
		G_replicated->g_length--;
		it++;
	}

	G_replicated->SetClaimVectors();
	G_replicated->pos_fork = 0;

	gnr_genes = G_replicated->gnr_genes;	//Copy gnr_genes to child before mutations have happened (and adjust during mutations).
	gnr_houses = G_replicated->gnr_houses;

	if (mutations_on)	//START mutations.
	{
		MutationList = new vector<bool>(g_length, false);	//Initialise the mutation-list with same length as the genome and with all entries set to 'false'.
		del_length = 0;
		dup_length = 0;
		pdup_length = &dup_length;
		pdel_length = &del_length;
		index = 0;

		//Mutate beads and set gnr_genes.
		it = BeadList->begin();
		while (it != BeadList->end())
		{
			if(IsGene(*it))
			{
				MutationList->at(index) = true;
			 	it = GeneMutate(it, pdel_length);
			}
			else if(IsTFBS(*it))
			{
				MutationList->at(index) = true;
				it = TFBSMutate(it, pdel_length);
			}
			else if(IsHouse(*it))
			{
				MutationList->at(index) = true;
				it = HouseMutate(it, pdel_length);
			}
			index++;
		}

		//Check that all beads had the chance to mutate.
		vector<bool>::iterator mit = MutationList->begin();
		while (mit != MutationList->end())
		{
			assert((*mit));
			mit++;
		}
		delete MutationList;
		MutationList = NULL;

		//Look for duplicated genes and tfbs's, which we will actually duplicate here.
		it = BeadList->begin();
		while (it != BeadList->end())
		{
			if ((*it)->type < 0)
			{
				if (IsGene(*it))
				{
					it=GeneDuplication(it, pdup_length);
				}
				else if (IsTFBS(*it))
				{
					it=TFBSDuplication(it);
					(*pdup_length)++;	//TFBS duplication always just adds one bead.
				}
				else if (IsHouse(*it))
				{
					it=HouseDuplication(it);
					(*pdup_length)++;
				}
			}
			else
			{
				it++;
			}
		}

		//Check that no more genes are tagged for duplication.
		it = BeadList->begin();
		while(it != BeadList->end())
		{
			assert ((*it)->type > 0);
			it++;
		}

		//Innovations of genes and TFBSs. Note that g_length is updated at GeneDuplication() and GeneDeletion() so any position along the genome is allowed for the novel gene/TFBS.
		if(uniform() < gene_innovation_mu)
		{
			it = GeneInnovation();
			PotentialTypeChange(it);
			(*pdup_length)++;
		}
		if(uniform() < tfbs_innovation_mu)
		{
			TFBSInnovation();
			(*pdup_length)++;
		}

		//Do shuffling mutations.
		it = BeadList->begin();
		while(it != BeadList->end())
		{
			if(IsGene(*it) && uniform() < gene_shuffle_mu)	it=GeneShuffle(it);
			else if(IsTFBS(*it) && uniform() < tfbs_shuffle_mu)	it=TFBSShuffle(it);
			else if(IsHouse(*it) && uniform() < house_shuffle_mu) it=HouseShuffle(it);
			else	it++;
		}

	}	//END of mutations.

	//Loop through the child genome and build the GeneTypes and GeneStates vectors.
	it = BeadList->begin();
	while (it != BeadList->end())
	{
		if(IsGene(*it))
		{
			gene = dynamic_cast<Gene*>(*it);
			index = FindIndexOfType(abs(gene->type));	//Find this type in GeneTypes vector.
			if(index == -1)	//Initialise this type in the GeneTypes and GeneStates vectors.
			{
				GeneTypes->push_back(abs(gene->type));
				GeneStates->push_back(gene->expression);
			}
			else
			{
				IncrementExpressionOfType(it);	//Find the index of this gene type in GeneTypes.
			}
		}
		it++;
	}

	if (g_length > 500)	//We only check after the full replication step, not after each replicated bead.
	{
		printf("Warning: genome sizes reached extravagant size (%d) during mitosis (i.e. due to mutations).\nExiting just to be safe...\n", g_length);
		cout << PrintContent(NULL, true, false) << endl;
		exit(1);
	}

	if(mutations_on)	assert(g_length == g_length_before_mut + (*pdup_length) - (*pdel_length));

	SetClaimVectors();	//First time the claim vectors of this child are set.

	pos_fork = 0;	//Set the fork to the beginning (nothing is replicated).
	pos_anti_ori = g_length;

}

void Genome::SplitGenome(Genome* G_replicated)	//Used to split a genome upon division
{
	//Find the fork with i_split.
	iter i_split = G_replicated->BeadList->begin();
	advance(i_split, G_replicated->pos_anti_ori);	//pos_anti_ori points to the end of the parental genome, which is now the first bead of the child genome.

	BeadList=new list<Bead*>();
	BeadList->splice(BeadList->begin(), *G_replicated->BeadList, i_split, G_replicated->BeadList->end());

	DevelopChildrenGenomes(G_replicated);
}


void Genome::ReplicateGenomeStep(int env, int res)
{
	iter it, start, end, it_new;
	Bead* bead;
	int gene_length = 0;
	int repl_remaining_steps = repl_step_size;

	repl_remaining_steps -= env;	//Subtract a number of steps defined by the current environment.

	if (repl_step_noise)
	{
		int repl_noise = ( (int)(uniform()*2) == 1 ) ? repl_remaining_steps:-repl_remaining_steps;
		repl_remaining_steps = ( (int)(uniform()*2) == 1 ) ? repl_remaining_steps:repl_remaining_steps+repl_noise;
	}

	repl_remaining_steps *= res;
	repl_remaining_steps /= 8;

	if (repl_remaining_steps <= 0)	return;

	start = BeadList->begin();
	advance(start, pos_fork);	//Now it points to the the first bead to be replicated in this replication step.
	end = start;	//End starts at start.

	//This loop sets "end".
	while ((distance(BeadList->begin(),end) < pos_anti_ori) && repl_remaining_steps > 0)	//The maximal position of end is defined by pos_anti_ori - 1 (pointing to the last bead of the parental genome). pos_anti_ori holds the number of genes in the parental genome, so if your distance to the first bead is pos_anti_ori, you are actually one past the last bead of the parental genome. The first replication step, this will point to NULL, but in consecutive steps it will point to a child bead; we want to point to an end point that does not change, hence the last bead of the parental genome.
	{
		if(IsGene(*end))	repl_remaining_steps--;	//Genes for sure count for the repl_step_size.
		else if(!replicate_entire_genes)	repl_remaining_steps--;	//TFBSs count if replicate_entire_genes is set to false.
		gene_length++;
		end++;
	}
	end--;

	bool last_round = false;
	it = start;	//Loop over a number of beads (how many we will replicate in one step).
	if(start == end)	last_round = true;	//We go straight into the last round.
	while (it != BeadList->end())
	{

		//Copy bead.
		bead=(*it)->Clone();
		(*BeadList).push_back(bead);
		g_length++;
		it_new = BeadList->end();
		it_new--;	//it_new now contains the just copied bead; needed to check whether it is a gene and do possible expression increment.
		if (IsGene(*it_new))
		{
			gnr_genes++;
			IncrementExpressionOfType(it_new);
		}
		else if (IsHouse(*it_new))
		{
			gnr_houses++;
		}

		it++;
		if (last_round)	break;
		if (it == end)	last_round = true;	//We have apparently hit the last bead of the parental genome, so time for one final replication step.
	}


	if (g_length > 500)	//We only check after the full replication step, not after each replicated bead.
	{
		printf("Warning: genome sizes reached extravagant size (%d) during replication.\nExiting just to be safe...\n", g_length);
		cout << PrintContent(NULL, true, false) << endl;
		exit(1);
	}

	pos_fork += gene_length;	//Move the fork to the right.
	if (pos_fork >= pos_anti_ori)	//If this is the final replication step.
	{
		pos_fork = pos_anti_ori;
		assert(g_length == 2*pos_anti_ori);
	}

	SetClaimVectors();	//Update ClaimVectors.
}

/*
###########################################################################
###########################################################################
							 |\ /\  |  | ----  /\  ---- o / \  |\ |
							|  V  \ l_J   |   /- \  |   | L_J  | \|
###########################################################################
###########################################################################
*/

Genome::iter Genome::GeneMutate(iter ii, int* pdel_len) {
	Gene* gene;
	gene = dynamic_cast<Gene*>(*ii);
	bool potential_type_change = false;
	int currentv, newv;

	double uu = uniform();
	if(uu < gene_duplication_mu)
	{
		(*ii)->type = -(*ii)->type;	//Mark for duplication during divison.
		mutant_genome = true;
		ii++;
	}

	else if(uu < gene_deletion_mu+gene_duplication_mu)
	{
		ii = GeneDeletion(ii, pdel_len);
		mutant_genome = true;
	}

	else
	{
		if(uniform() < gene_threshold_mu)	//All mutational events are now independent, i.e. a gene can get multiple mutations or none, or one.
		{
			if (uniform()>0.8)	gene->threshold = (uniform()>0.5) ? gene->threshold+1 : gene->threshold-1;
			else
			{
				currentv = gene->threshold;
				newv = gene->threshold;
				while(newv == currentv)
				{
					gene->threshold = (int)(uniform()*(2*WeightRange+1) - WeightRange);
					newv = gene->threshold;
				}
			}
			mutant_genome = true;
		}

		if (uniform() < gene_activity_mu)
		{
			if (uniform()>0.8)	gene->activity = (uniform()>0.5) ? gene->activity+1 : gene->activity-1;
			else
			{
				currentv = gene->activity;
				newv = gene->activity;
				while(newv == currentv)
				{
					gene->activity = (uniform()>0.5) ? -1*uniform()*WeightRange : 1*uniform()*WeightRange;
					newv = gene->activity;
				}
			}
			potential_type_change = true;
			mutant_genome = true;
		}

		for(int k=0; k<binding_length; k++)
		{
			if(uniform() < gene_binding_domain_mu)
			{
				if (gene->binding_domain[k] == false) gene->binding_domain[k] = true;
				else if (gene->binding_domain[k] == true) gene->binding_domain[k] = false;
			 	potential_type_change = true;
				mutant_genome = true;
			}
		}

		if (potential_type_change){
			PotentialTypeChange(ii);	//Check that it has not become the same type as one of the other genes. I don't think it matters that it is one of the original five. Every gene can convert into an existing type; for duplicated genes, one can attain a new type.
		}
		ii++;
	}

	return ii;
}

Genome::iter Genome::TFBSMutate(iter ii, int* pdel_len)
{
	TFBS* tfbs;
	tfbs = dynamic_cast<TFBS*>(*ii);		//downcast the Bead-pointer to a TFBS-pointer with the name tfbs;
	int currentv, newv;

	double uu = uniform();
	if(uu < tfbs_duplication_mu)
	{
		(*ii)->type = -(*ii)->type;
		mutant_genome = true;
		ii++;
	}

	else if(uu < tfbs_duplication_mu+tfbs_deletion_mu)
	{
		ii = TFBSDeletion(ii);
		(*pdel_len)++;
		mutant_genome = true;
	}

	else
	{
		for (int k=0; k<binding_length; k++)
		{
			if(uniform() < tfbs_binding_site_mu)
			{
				if (tfbs->binding_site[k] == false) tfbs->binding_site[k] = true;
				else if (tfbs->binding_site[k] == true) tfbs->binding_site[k] = false;
				mutant_genome = true;
			}
		}
		if(uniform() < tfbs_activity_mu)
		{
			if (uniform()>0.8)	tfbs->activity = (uniform()>0.5) ? tfbs->activity+1 : tfbs->activity-1;
			else
			{
				currentv = tfbs->activity;
				newv = tfbs->activity;
				while(newv == currentv)
				{
					tfbs->activity = (uniform()>0.5) ? -1*uniform()*WeightRange : 1*uniform()*WeightRange;
					newv = tfbs->activity;
				}
			}
			mutant_genome = true;
		}
		ii++;
	}
	return ii;
}

Genome::iter Genome::HouseMutate(iter ii, int* pdel_len)
{
	double uu = uniform();

	if(uu < house_duplication_mu)
	{
		// cout << "Hous dup" << endl;
		(*ii)->type = -(*ii)->type;
		mutant_genome = true;
		ii++;
		// cout << "House duped" << endl;
	}
	else if(uu < house_duplication_mu+house_deletion_mu)
	{
		// cout << "HOuse del" << endl;
		ii = HouseDeletion(ii);
		(*pdel_len)++;
		mutant_genome = true;
	}
	else	ii++;

	return ii;
}

Genome::iter Genome::GeneDuplication(iter ii, int* pdup_len)
{
	int copy_length;
	iter insertsite, first, last;
	list<Bead*> BeadListTemp;	//Create a new temporary genome list.

	(*ii)->type = -(*ii)->type;	//Remove the duplication flag.

	//Copy the gene with its upstream tfbs's to a temporary chromosome.
	last = ii;
	last++;   //One further than the gene position (the one not to be part of the dupl).
	first=FindFirstTFBSInFrontOfGene(ii);	//First tfbs in front of gene.
	copy_length=distance(first, last);
	CopyPartOfGenomeToTemplate(first, last, &BeadListTemp); //Makes a 'chromosome' with only this gene (and its tfbs) on it.

	//Splice the temporary chromosome into the full genome.
	insertsite=FindRandomGenePosition();			//Find position of gene to insert in front of.
	insertsite=FindFirstTFBSInFrontOfGene(insertsite);	//Find first tfbs in front of this gene.
	BeadList->splice(insertsite, BeadListTemp);	//Splice temporary list into chromosome.

	//Increment the number of beads and the number of genes.
	g_length+=copy_length;
	(*pdup_len)+=copy_length;
	gnr_genes++;

	ii=last;	//Make sure ii points to one position further than just duplicated gene.
	return ii;
}

Genome::iter Genome::GeneInnovation()
{
	Gene* gene;
	iter insertsite;

	gene = new Gene();
	gene->RandomGene();

	insertsite=FindRandomGenePosition();
	insertsite=FindFirstTFBSInFrontOfGene(insertsite);
	insertsite=BeadList->insert(insertsite, gene);

	g_length++;
	gnr_genes++;
	return insertsite;
}

Genome::iter Genome::GeneDeletion(iter ii, int* pdel_len)
{
	iter first, last, jj;
	int del_length;

	last=ii;//gene position
	last++;//one further than the gene position
	first=FindFirstTFBSInFrontOfGene(ii);//first tfbs in front of gene

	//Decrement the number of beads and the number of genes.
	del_length = distance(first, last);
	g_length -= del_length;
	gnr_genes--;
	(*pdel_len) += del_length;

	jj=first;
	while( jj != last )
	{
		delete *jj;		//Here you delete the element that the iterator is pointing to.
		jj++;
	}
	ii=(*BeadList).erase(first, last);		//Here you remove the iterators, i.e. the pointers to the deleted elements from the iterator 'list'.
	return ii;
}

Genome::iter Genome::GeneShuffle(iter ii)
{
	iter insertsite, first, last, jj;
	list<Bead*> BeadListTemp;	//Create a new temporary genome list.

	//Copy the gene with its upstream tfbs's to a temporary chromosome.
	last = ii;
	last++;   //One further than the gene position (the one not to be part of the dupl).
	first=FindFirstTFBSInFrontOfGene(ii);	//First tfbs in front of gene.
	CopyPartOfGenomeToTemplate(first, last, &BeadListTemp); //Makes a 'chromosome' with only this gene (and its tfbs) on it.
	last--;	//This is important, because if we move the virtual copy of the gene directly downstream of its original (i.e. it does not really move), then if last is still pointing to the bead that was originally adjacent to the gene, both the virtual copy and the original will be removed.

	//Splice the temporary chromosome into the full genome.
	insertsite=FindRandomGenePosition();			//Find position of gene to insert in front of.
	insertsite=FindFirstTFBSInFrontOfGene(insertsite);	//Find first tfbs in front of this gene.
	BeadList->splice(insertsite, BeadListTemp);	//Splice temporary list into chromosome.

	last++;	//Now make last point to the next bead after the original gene again (see above), so that we will only remove the original gene with its preceding binding sites.
	//Remove the bead from its original position, taken from GeneDeletion().
	jj=first;
	while( jj != last )
	{
		delete *jj;
		jj++;
	}
	ii=(*BeadList).erase(first, last);
	return ii;
}

Genome::iter Genome::TFBSDuplication(iter ii)
{
	iter tt, upstream;
	int randpos;
	TFBS* tfbs;
	tfbs=dynamic_cast<TFBS *>(*ii);

	tfbs->type = -tfbs->type;	//Make the type positive, so that we know that the just duplicated TFBS has been randomly moved to a new position on the genome.
	TFBS* tfbsnew = new TFBS(*tfbs);

	tt = (*BeadList).begin();
	randpos = (int)(uniform()*g_length);
	advance(tt,randpos);			// tt holds random spot in the genome e.g. |---------------x-----------|
	tt = (*BeadList).insert(tt, tfbsnew);	//Insert tfbs-copy to the left of a random position in the genome (tt).

	g_length++;
	ii++;
	return ii;
}

void Genome::TFBSInnovation()
{
	TFBS* tfbs;
	iter insertsite;
	int randpos;

	tfbs = new TFBS();
	tfbs->RandomTFBS();

	randpos = (int)(uniform() * g_length);
	insertsite = BeadList->begin();
	advance(insertsite,randpos);
	BeadList->insert(insertsite, tfbs);

	g_length++;
}

Genome::iter Genome::TFBSDeletion(iter ii)
{
	TFBS *tfbs;
	tfbs=dynamic_cast<TFBS *>(*ii);
	delete (tfbs);
	ii=(*BeadList).erase(ii);

	g_length--;
	return ii;
}

Genome::iter Genome::TFBSShuffle(iter ii)
{
	iter tt, upstream;
	int randpos;
	TFBS* tfbs;
	tfbs=dynamic_cast<TFBS *>(*ii);

	//Create copy and insert at random location.
	TFBS* tfbsnew = new TFBS(*tfbs);

	tt = (*BeadList).begin();
	randpos = (int)(uniform()*g_length);
	advance(tt,randpos);			// tt holds random spot in the genome e.g. |---------------x-----------|
	tt = (*BeadList).insert(tt, tfbsnew);	//Insert tfbs-copy to the left of a random position in the genome (tt).

	//Remove the old bead.
	delete (tfbs);
	ii=(*BeadList).erase(ii);

	return ii;
}

Genome::iter Genome::HouseDuplication(iter ii)
{
	iter tt, upstream;
	int randpos;

	(*ii)->type = -(*ii)->type;	//Make the type positive, so that we know that the just duplicated TFBS has been randomly moved to a new position on the genome.
	House* housenew = new House();

	tt = (*BeadList).begin();
	randpos = (int)(uniform()*g_length);
	advance(tt,randpos);			// tt holds random spot in the genome e.g. |---------------x-----------|
	tt = (*BeadList).insert(tt, housenew);	//Insert tfbs-copy to the left of a random position in the genome (tt).

	gnr_houses++;
	g_length++;
	ii++;
	return ii;
}

Genome::iter Genome::HouseDeletion(iter ii)
{
	delete (*ii);
	ii=(*BeadList).erase(ii);

	gnr_houses--;
	g_length--;
	return ii;
}

Genome::iter Genome::HouseShuffle(iter ii)
{
	iter tt, upstream;
	int randpos;
	House* house;
	house = dynamic_cast<House *>(*ii);

	//Create copy and insert at random location.
	House* housenew = new House();

	tt = (*BeadList).begin();
	randpos = (int)(uniform()*g_length);
	advance(tt,randpos);			// tt holds random spot in the genome e.g. |---------------x-----------|
	tt = (*BeadList).insert(tt, housenew);	//Insert tfbs-copy to the left of a random position in the genome (tt).

	//Remove the old bead.
	delete (house);
	ii=(*BeadList).erase(ii);

	return ii;
}

Genome::iter Genome::FindFirstTFBSInFrontOfGene(iter ii) const
{
	reviter rii(ii);	//BEHIND in this case means an element to the left in the list (i.e if ii points to the 6th element, rii(ii) will make rii point to the 5th element). The function .base() as used below will make rii.base() point to the 6th element again.

	reviter jj = (*BeadList).rend();//search should be bounded
	while(rii != jj)//begin not yet reached
	{
		if(IsTFBS(*rii))	rii++;
		else	jj = rii;
	}
	return jj.base();
}

Genome::iter Genome::FindRandomGenePosition() const
{
	std::list< iter > pos;
	std::list< iter >::iterator ipos;
	iter i, ii;
	int randpos;

	if(gnr_genes+gnr_houses==0)	return (*BeadList).end();
	else
	{
		i=(*BeadList).begin();
		while(i != (*BeadList).end())
		{
			if(IsGene(*i) || IsHouse(*i))	pos.push_back(i);
			i++;
		}
		randpos=(int)(uniform()*(gnr_genes+gnr_houses));	//If you found the first gene, randpos will be 0 (no need to advance to another gene); if you find the last gene, randpos will be gnr_genes-1 which is enough to reach the gnr_genes'th gene.
		ii = (*boost::next(pos.begin(),randpos));	//Possibly try advance(ii, randpos) instead.
		return ii;
	}
}

void Genome::PotentialTypeChange(iter ii)
{
	bool found_matching_type = false;
	int type_abundance = CountTypeAbundance(abs((*ii)->type));
	list<int> UsedTypes;

	iter jj = BeadList->begin();	//The other genes in the genome
	while(jj != BeadList->end())
	{
		if(IsGene(*jj) && jj!=ii)	//We don't convert genes to themselves.
		{
			UsedTypes.push_back(abs((*jj)->type));
			bool genes_are_the_same = CheckSameGeneTypes(ii, jj);

			if(genes_are_the_same)
			{
				(*ii)->type = abs((*jj)->type);		//Convert to existing gene type.
				// IncrementExpressionOfType(ii);	//We could here increment the expression of the new type if it was expressed or something
				found_matching_type = true;
				return;		//We have found a match, converted the gene; time to try mutation of the next bead.
			}
		}
		jj++;
	}
	if(found_matching_type == false && (type_abundance > 1 || (*ii)->type==0))	//We haven't been able to convert it to an existing type, so let's define it as a new type. type_abundance should always be more than 1 because there is always an original copy on the parental section of the genome.
	{
		list<int>::iterator tit;

		for(int x=1; x<1001; x++)		//Find a not yet used number to use as the type.
		{
			//Again go through all beads to find an unused type.
			tit = find(UsedTypes.begin(), UsedTypes.end(), x);
			if(tit == UsedTypes.end())	//X is not yet used as a gene type.
			{
				(*ii)->type = x;
				return;
			}
		}
	}
}

bool Genome::CheckSameGeneTypes(iter ii, iter jj)
{
	Gene* gene_ii = dynamic_cast<Gene*>(*ii);
	Gene* gene_jj = dynamic_cast<Gene*>(*jj);

	bool genes_are_the_same = true;		//Prove me wrong!
	if(gene_ii->activity != gene_jj->activity) genes_are_the_same = false;
	for (int k=0; k<binding_length; k++)
	{
		if(gene_ii->binding_domain[k] != gene_jj->binding_domain[k]) genes_are_the_same = false;
	}
	return genes_are_the_same;
}

bool Genome::CheckSameGeneTypes(Gene* gene_ii, Gene* gene_jj)
{
	bool genes_are_the_same = true;		//Prove me wrong!
	if(gene_ii->activity != gene_jj->activity) genes_are_the_same = false;
	for (int k=0; k<binding_length; k++)
	{
		if(gene_ii->binding_domain[k] != gene_jj->binding_domain[k]) genes_are_the_same = false;
	}
	return genes_are_the_same;
}

int Genome::CountTypeAbundance(int type)
{
	int count_type = 0;
	iter ii = BeadList->begin();
	while (ii != BeadList->end())
	{
		if((*ii)!=NULL && IsGene(*ii))
		{
			if(abs((*ii)->type)==type)	count_type++;
		}
		ii++;
	}
	return count_type;
}

/*
###########################################################################
###########################################################################
(((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))
(((((((((((((((((((((((((((((((((((((())))))))))))))))))))))))))))))))))))))
###########################################################################
###########################################################################
*/

void Genome::ReadInitialGenome()
{
	BeadList = new list<Bead*>();
	GeneTypes = new vector<int>();
	GeneStates = new vector<int>();	//If we set -g [G1,S,G2,M] in the command line, we will update genestates while we read the genome.

	ifstream infile(genome_init.c_str());
	string line;

	if (!infile.is_open())
	{
		printf("Genome input file could not be opened.\n");
		exit(1);
	}

	printf("Reading genome from file: %s\n", genome_init.c_str());
	while(getline(infile,line))
	{
		ReadBeadsFromString(line);
	}

	//Read or randomly set initial expression.
	if(genestate_init == "")	for(int g=0; (size_t)g<GeneTypes->size(); g++)	GeneStates->push_back((int)(uniform()*2));
	else if(genestate_init != "G1" && genestate_init != "S" && genestate_init != "G2" && genestate_init != "M")	ReadInitialGeneStates();	//if genestate_input was specified, the GeneStates should have already been initialised during the reading of the genome.

	//Set forks.
	pos_fork = 0;
	pos_anti_ori = g_length;

	SetClaimVectors();
}


void Genome::ReadBeadsFromString(string genome)
{
	char* bead;
	int index;
	bool StageInit[5];
	gnr_genes = 0;
	gnr_houses = 0;
	g_length = 0;
	Gene* gene;
	TFBS* tfbs;
	House* house;
	int type, threshold, activity, q;
	char* buffer;
	bool bitstring[binding_length];

	if (genestate_init == "G1")	for(int g=1; g<6; g++)	StageInit[g] = StageTargets[0][g-1];
	else if (genestate_init == "S") for(int g=1; g<6; g++)	StageInit[g] = StageTargets[1][g-1];
	else if (genestate_init == "G2") for(int g=1; g<6; g++)	StageInit[g] = StageTargets[2][g-1];
	else if (genestate_init == "M") for(int g=1; g<6; g++)	StageInit[g] = StageTargets[3][g-1];
	else if (backup_reboot == "")	printf("Expression not set to a cell-cycle stage.\n");

	bead = strtok((char*)genome.c_str(),".");
	while (bead != NULL)
	{
		if(bead[1] == 'G')	//Bead is a gene
		{
			buffer = new char[binding_length+2];	//The character array needs to be longer than the actual number of characters is should hold. In this case (and the one below) I think it also happens that the last ')' is stored in the array; which is nice because I use it below to know that bitstring is finished.
			int success = sscanf(bead, "(G%d:%d:%d:%s)", &type, &threshold, &activity, buffer);
			if(success != 4) cerr << "Could not find sufficient information for this gene. Genome file potentially corrupt. \n" << endl;
			q = 0;
			while(buffer[q] != ')')		//The extra bracket stored in buffer actually pays off here because I don't know how else I would know that we reached the end of buffer.
			{
				bitstring[q] = (buffer[q]=='1');	//Easiest way I could think of to convert a character to a boolean; (bool)buffer[q] always returns 1 (whether '1' or '0')!
				q++;
			}
			gene = new Gene(type, threshold, activity, bitstring, (abs(type)<6));	//We initialise with zero expression, but these should be updated with the first round of UpdateGeneStates().
			index = FindIndexOfType(abs(type));
			if(index == -1)	//Type not found in GeneTypes, so add an element.
			{
				GeneTypes->push_back(abs(type));
				if (genestate_init == "G1" || genestate_init == "S" || genestate_init == "G2" || genestate_init == "M")	//Then we have given a cell-cycle stage as input argument so let's set the five main TF's to their corresponding value.
				{
					if (abs(type) < 6)	//For the five main TFs we set the expression according to the stage given as input argument.
					{
						if(StageInit[abs(type)])	GeneStates->push_back(1);
						else	GeneStates->push_back(0);
					}
					else
					{
						GeneStates->push_back(0);
						// GeneStates->push_back((int)(uniform()*2));	//For the genes other than the five main TFs we randomly set the expression.
					}
				}
			}
			else if (genestate_init == "G1" || genestate_init == "S" || genestate_init == "G2" || genestate_init == "M"){
				//This is the second copy of a particular type, so just increment the existing GeneStates element.
				//This gene type is expressed if it is active in the specified cell-cycle stage.
				if (StageInit[abs(type)])	GeneStates->at(index) += 1;
				else if(abs(type) > 6)	GeneStates->at(index) += (int)(uniform()*2);	//Increment gene state by 0 or 1.
			}
			(*BeadList).push_back(gene);
			gnr_genes++;
			g_length++;
			delete [] buffer;	//I think this prevents a memory leak. Otherwise buffer simply frees up new memory the next time it says "buffer = new char()". Because I did not use "new []" to create the variable, I should not use "delete []" to free up the memory.
			buffer = NULL;
		}
		else if (bead[1] == 'H')	//Bead is a house
		{
			house = new House();
			(*BeadList).push_back(house);
			gnr_houses++;
			g_length++;
		}
		else	//Bead is a tfbs
		{
			buffer = new char[binding_length+2];
			int success = sscanf(bead, "(%d:%s)", &activity, buffer);
			if(success != 2) cerr << "Could not find sufficient information for this TFBS. Genome file potentially corrupt. \n" << endl;
			q = 0;
			while(buffer[q] != ')')
			{
				bitstring[q] = (buffer[q]=='1');
				q++;
			}
			tfbs = new TFBS(1, activity, bitstring);
			(*BeadList).push_back(tfbs);
			g_length++;
			delete [] buffer;
			buffer = NULL;
		}
		bead = strtok(NULL, ".");
	}
}


void Genome::ReadInitialGeneStates()
{
	ifstream infile(genestate_init.c_str());
	string line;
	char* token;

	if (!infile.is_open())
	{
		printf("Gene-state input file could not be opened.\n");
		exit(1);
	}

	printf("Reading gene states from file: %s\n", genestate_init.c_str());
	while(getline(infile,line))
	{
		string subline = line.substr(1,line.size()-2);
		token = strtok((char*)subline.c_str(),", ");
		for(int g=0; (size_t)g<GeneTypes->size(); g++)
		{
			if (token == NULL)
			{
				cerr << "Error: Gene state file did not contain enough states." << endl;
				exit(1);
			}
			GeneStates->push_back(atoi(token));
			token = strtok(NULL, ", ");
		}
	}
}

void Genome::InitialiseRandomGenome()
{
	BeadList = new list<Bead*>();	//Create the empty beadlist for this genome.
	gnr_genes = 0;
	gnr_houses = 0;
	g_length = 0;

	Gene* gene;
	TFBS* tfbs;
	House* house;

	//Just place all household genes at the beginning for now.
	for (int houses=1; houses<=nr_household_genes; houses++)
	{
		house = new House();
		(*BeadList).push_back(house);
		gnr_houses++;
		g_length++;
	}

	//For now, I don't mind that genomes are ordered by gene type.
	for (int gene_type=1; gene_type<=init_nr_gene_types; gene_type++)
	{

		//Define its tfbs's.
		for (int j=0; j<init_nr_tfbs_per_gene; j++)
		{
			// activity = (uniform()>0.5) ? -1 : 1;
			// for (int k=0; k<binding_length; k++)	binding_section[k] = (uniform()>0.5) ? true : false;
			// tfbs = new TFBS(tfbs_type, activity, binding_section);
			tfbs = new TFBS();
			tfbs->RandomTFBS();
			(*BeadList).push_back(tfbs);
			g_length++;
		}

		//Define the gene itself. Here the type matters, and this is the only case where we first want to say that it is a particular type and then make sure that it has a unique bitstring. This bitstring may mutate, in which case we should track it to know whether it is one of the five main types. Conversely, I could say that the five main types are always defined by five particular genes. If these genes the definition of the type also changes. Other genes will be of different types until they happen to randomly attain the same properties as one of the original five genes, in which case they are set to that type. All this is important for duplication and deletion events.
		// threshold = (int)(uniform()*(2*WeightRange+1) - (int)WeightRange);	//Threshold between -3 and 3 (including these borders).
		// activity = (uniform()>0.5) ? -1 : 1;
		// for (int k=0; k<binding_length; k++)	binding_section[k] = (uniform()>0.5) ? true : false;
		gene = new Gene();
		gene->RandomGene();

		//Check that we have not initiated a gene of the same type
		iter jj = BeadList->begin();
		while(jj != BeadList->end())
		{
			if(IsGene(*jj))
			{
				Gene* check_gene = dynamic_cast<Gene*>(*jj);
				bool genes_are_the_same = CheckSameGeneTypes(gene, check_gene);

				if(genes_are_the_same)
				{
					for(int t=0; t<init_nr_tfbs_per_gene; t++)	(*BeadList).pop_back();
					// g_length -= init_nr_gene_types;
					g_length -= init_nr_tfbs_per_gene;
					delete gene;
					gene_type--;
					break;
				}
			}
			jj++;
		}

		gene->type = gene_type;
		// gene = new Gene(gene_type, threshold, activity, binding_section, 0);	//We should only get here if the new gene is not the same as one of the old genes.
		(*BeadList).push_back(gene);
		gnr_genes++;
		g_length++;

	}
	//Randomly initialise gene expression [0,g] for each type (where g is the number of genes of that type).
	GeneStates = new vector<int>();
	GeneTypes = new vector<int>();
	for(int g=1; g<=init_nr_gene_types; g++)
	{
		GeneStates->push_back((int)(uniform()*2));
		GeneTypes->push_back(g);
	}

	pos_fork = 0;
	pos_anti_ori = g_length;

	SetClaimVectors();	//For each TFBS, put all bindings strengths with all genes in its claim vector.
}

int Genome::MatchNextState(int state)
{
	int expression, index, match_next_state = 0;
	Genome::gene_iter git;

	for (int g=1; g<6; g++)	//g is the type of gene that we are looking for, ie G1-G5; but these may be shuffled (or some missing) from the actual GeneStates.
	{
		git = find(GeneTypes->begin(), GeneTypes->end(), g);
		if (git == GeneTypes->end())	expression = 0;
		else
		{
			index = distance(GeneTypes->begin(), git);
			expression = GeneStates->at(index);
		}

		if(  (expression==0 && !StageTargets[state][g-1])  ||  (expression!=0 && StageTargets[state][g-1])  )
		{
			match_next_state++;
		}
		// else	break;	//If one gene is not expressed correctly, there is no use in looking at the other genes.
	}
	return match_next_state;
}

void Genome::UpdateGeneStates()
{
	TFBS* tfbs;
	Gene* gene;

	vector<int>* NextStateExpression;
	NextStateExpression = new vector<int>();
	NextStateExpression->resize(GeneStates->size());

	for(int g=0; (size_t)g<GeneStates->size(); g++)	NextStateExpression->at(g) = 0;	//Make sure vector is empty (i.e. all zeroes).
	int sumeffects=0;
	iter it = BeadList->begin();
	while (it != BeadList->end())
	{
		if(IsGene(*it))
		{
			gene = dynamic_cast<Gene*>(*it);

			sumeffects -= gene->threshold;
			int index = FindIndexOfType(abs(gene->type));	//Find the right index for GeneTypes/GeneStates of the gene type that we have in our hands.
			gene->expression = max(min(sumeffects+1,1),0);	//For the +1, see explanation of the gene threshold in Gene.hh
			NextStateExpression->at(index) += gene->expression;
			sumeffects = 0;
		}
		else if(IsTFBS(*it))
		{
			iter i_gene;
			i_gene = MatchGeneToTFBS(it);		//We actually choose a specific gene of the winning gene type / or we pick one randomly since they will all give the same activity.
			if (i_gene == BeadList->end())
			{
				sumeffects += 0;
			}
			else
			{
				gene = dynamic_cast<Gene*>(*i_gene);
				tfbs = dynamic_cast<TFBS*>(*it);
				sumeffects += tfbs->activity*gene->activity;
			}
		}
		it++;
	}

	//Now set the NextStateExpression into the GeneStates vector.
	gene_iter inext = NextStateExpression->begin();
	gene_iter irenew = GeneStates->begin();
	while (inext != NextStateExpression->end())
	{
		(*irenew) = (*inext);	//Put the value of the NextStateExpression in the GeneStates.
		inext++;
		irenew++;
	}
	//Now release the memory from NextStateExpression.
	delete NextStateExpression;
	NextStateExpression = NULL;	//Not sure if this is necessary.
}

void Genome::SetClaimVectors()
{
	iter it, ig;
	it = BeadList->begin();
	while(it != BeadList->end())
	{
		if(IsTFBS(*it))
		{
			//Return its binding strength for every gene in the genome.
			TFBS* tfbs;
			tfbs = dynamic_cast<TFBS*>(*it);
			tfbs->ClaimVector->resize(GeneStates->size());

			ig = BeadList->begin();
			while(ig != BeadList->end())
			{
				if(IsGene(*ig))		//Currently we are still recalculating the bitstring match for each gene in our genome, even though genes of the same type should only have to be matched once.
				{
					int index = FindIndexOfType(abs((*ig)->type));
					tfbs->ClaimVector->at(index) = '0'+MatchBitStrings(*it, *ig);
				}
				ig++;
			}
		}
		it++;
	}
}


Genome::iter Genome::MatchGeneToTFBS(iter i_tfbs)
{
	TFBS* tfbs;
	tfbs = dynamic_cast<TFBS*>(*i_tfbs);
	double claim_value, final_claim;
	double claim_sum = 0.0;
	iter i_gene;
	//Calculate the claim sum.	//Maybe faster to iterate through the elements.
	for(int g=0; (size_t)g<GeneStates->size(); g++)
	{
		claim_value = (int)tfbs->ClaimVector->at(g) - 48;
		claim_sum += pow( ((double)(claim_value) / (double)(binding_length)), tfbs_selection_exponent) * (double) GeneStates->at(g);	//Claim of gene types is proportional to the number of active genes of that type and to their binding strength to the tfbs.
	}
	//Based on the cumulative claim, pick one of the genes.
	i_gene = BeadList->begin();
	if (claim_sum == 0.0)  return BeadList->end();	//There is no active TF that can bind the TFBS.
  else
	{
		double total_claim = max(empty_tf_claim_zero,claim_sum);
		double die_roll = uniform()*total_claim;

		double empty_claim = max(0., empty_tf_claim_zero-claim_sum);	//Do this to quickly finish TFBS-TF matching in case no active TF managed to bind.
		if(die_roll <= empty_claim)	return BeadList->end();	//No active TF managed to bind the TFBS.
		else	die_roll -= empty_claim;
		for(int g=0; (size_t)g<GeneStates->size(); g++)
		{
			claim_value = (int)tfbs->ClaimVector->at(g) - 48;
			final_claim = pow( ((double)(claim_value) / (double)(binding_length)), tfbs_selection_exponent) * (double)(GeneStates->at(g));
			if ( die_roll <= final_claim )
			{
				//Increment i_gene until we find a gene of the correct type.
				while(!(IsGene(*i_gene) && abs((*i_gene)->type)==GeneTypes->at(g)))
				{
					i_gene++;
				}
				return i_gene;
			}
			die_roll -= final_claim;
		}
	}
	printf("Error: claim out of bounds.\n");
	exit(1);
}

int Genome::MatchBitStrings(Bead* b_tfbs, Bead* b_gene){
	TFBS* tfbs;
	Gene* gene;
	tfbs = dynamic_cast<TFBS*> (b_tfbs);
	gene = dynamic_cast<Gene*> (b_gene);

  int binding_strength = 0;
  for(int k=0; k<binding_length; k++)
	{
    if(tfbs->binding_site[k] != gene->binding_domain[k]) binding_strength++;		//We do complementary pairing for the beauty of it (but also expected in visualisations).
  }
	return binding_strength;
}

/*
###########################################################################
###########################################################################
 |\  /\  O  r-- /--  /_-  |   |
|  V  |  |  _\  \__  \_  |__ |__
###########################################################################
###########################################################################
*/

void Genome::GenomeToNetwork(double** Network)
{
	double NetRow[8] = {.0, .0, .0, .0, .0, .0, .0, .0};	//First element is threshold of gene (sum of the copies), then 5 elements for incoming effects of G1-5, and then 2 elements for other positive incoming effects (sum) and other negative incoming effects (sum).
	vector<int> TypeCopyNumber;
	vector<int> GeneActivity;
	iter it;
	gene_iter git;
	int index, type, gene_act, gene_copynum, bead_count = 0, claim_value;
	TFBS* tfbs;
	TFBS::claim_iter cit;
	double Effect;
	Gene* gene;

	//Before we can interpret the ClaimVectors of all TFBS's we need to know if there are genes present in higher copy number.
	TypeCopyNumber.resize(GeneTypes->size(),0);
	GeneActivity.resize(GeneTypes->size(),0);

	it = BeadList->begin();
	while(it != BeadList->end() && bead_count < pos_anti_ori)
	{
		if(IsGene(*it))
		{
			git = find(GeneTypes->begin(), GeneTypes->end(), (*it)->type);
			index = distance(GeneTypes->begin(), git);
			TypeCopyNumber.at(index)++;
			gene = dynamic_cast<Gene*>(*it);
			GeneActivity.at(index) = gene->activity;
		}
		it++;
		bead_count++;
	}

	bead_count = 0;
	it = BeadList->begin();
	while(it != BeadList->end() && bead_count < pos_anti_ori)
	{
		if(IsTFBS(*it))
		{
			//Iterate through ClaimVector, noting whether it is a claim of G1-5 or a different gene.
			tfbs = dynamic_cast<TFBS*>(*it);
			cit = tfbs->ClaimVector->begin();
			index = 0;
			while (cit != tfbs->ClaimVector->end())
			{
				type = GeneTypes->at(index);
				gene_act = GeneActivity.at(index);
				gene_copynum = TypeCopyNumber.at(index);

				claim_value = (int)(*cit) - 48;
				Effect = pow( ((double)(claim_value) / (double)(binding_length)), tfbs_selection_exponent) * gene_copynum * gene_act * tfbs->activity;	//(*cit) is the raw claim, i.e. the bitstring match between TFBS and gene.

				if(type <= 5)	NetRow[type] += Effect;
				else if (Effect > 0)	NetRow[6] += Effect;
				else	NetRow[7] += Effect;

				cit++;
				index++;
			}
		}
		else if(IsGene(*it))
		{
			if((*it)->type <= 5)
			{
				gene = dynamic_cast<Gene*>(*it);
				NetRow[0] += gene->threshold;

				//Put NetRow in the right row of Net.
				for (int col=0; col<8; col++)
				{
					Network[gene->type-1][col] += NetRow[col];
					NetRow[col] = 0.;
				}
			}
		}
		it++;
		bead_count++;
	}
}


bool Genome::IsGene(Bead* bead) const	{
	return (bool)(typeid(*bead) == typeid(Gene));	// typeid determines the class of an object at runtime
}

bool Genome::IsTFBS(Bead* bead) const {
	return (bool)(typeid(*bead) == typeid(TFBS));
}

bool Genome::IsHouse(Bead* bead) const {
	return (bool)(typeid(*bead) == typeid(House));
}

void Genome::DecrementExpressionOfType(iter ii)
{
	Gene* gene = dynamic_cast<Gene*>(*ii);
	int index = FindIndexOfType(abs(gene->type));	//abs() makes sure that genes flagged for duplication by a '-' are not treated as different types.

	GeneStates->at(index) -= gene->expression;
}

void Genome::IncrementExpressionOfType(iter ii)
{
	Gene* gene = dynamic_cast<Gene*>(*ii);
	int index = FindIndexOfType(abs(gene->type));

	GeneStates->at(index) += gene->expression;
}

int Genome::FindIndexOfType(int type)
{
	gene_iter git = find(GeneTypes->begin(), GeneTypes->end(), type);
	if(git == GeneTypes->end())	return -1;	//This is a signal during MoveGenomeToChild() that the type does not yet exist in the GeneTypes vector.
	else	return distance(GeneTypes->begin(), git);
}

string Genome::PrintContent(list<Bead*> * chromosome, bool terminal, bool only_parent)	// Printing function for terminal output.
{
	string GenomeContent;
	if(chromosome == NULL) chromosome = this->BeadList;
	iter i = chromosome->begin();
	Gene *gene;
	TFBS *tfbs;

	iter end;
	if (only_parent)
	{
		end = chromosome->begin();
		advance(end, pos_anti_ori);
	}
	else	end = chromosome->end();

	string gene_color_prefix = "";
	string gene_color_suffix = "";
	string tfbs_color_prefix = "";
	string tfbs_color_suffix = "";
	string house_color_prefix = "";
	string house_color_suffix = "";

	if(terminal==true){
		gene_color_prefix = "\033[94m";
		gene_color_suffix = "\033[0m";
		tfbs_color_prefix = "\033[92m";
		tfbs_color_suffix = "\033[0m";
		house_color_prefix = "\033[95m";
		house_color_suffix = "\033[0m";
	}

	while(i!=end) {
		std::stringstream stringtemp;
		if(i!=chromosome->begin()) GenomeContent +=".";
		GenomeContent += "(";

		if(IsGene(*i)) {
			gene=dynamic_cast<Gene *>(*i);
			stringtemp << gene_color_prefix << "G" << gene->type << ":" << gene->threshold << ":" << gene->activity << ":";
			for(int k=0; k<binding_length; k++)	stringtemp << gene->binding_domain[k];
			stringtemp << gene_color_suffix;
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsTFBS(*i)) {
			tfbs=dynamic_cast<TFBS *>(*i);
			stringtemp << tfbs_color_prefix << tfbs->activity << ":";
			for(int k=0; k<binding_length; k++)	stringtemp << tfbs->binding_site[k];
			stringtemp << tfbs_color_suffix;
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}
		else if(IsHouse(*i)) {
			stringtemp << house_color_prefix << "H" << house_color_suffix;
			GenomeContent+=stringtemp.str();
			stringtemp.clear();
		}

		GenomeContent += ")";
		i++;
	}
	return GenomeContent;
}

string Genome::PrintGeneStateContent(bool sorted)
{
	gene_iter git;
	int index, count_types=0;
	string GeneStateContent = "[";

	if (!sorted)
	{
		for (int i=0; (size_t)i<GeneStates->size(); i++)
		{
			stringstream stringtemp;
			stringtemp << GeneStates->at(i);
			if((size_t)i<GeneStates->size()-1) GeneStateContent += stringtemp.str() + ", ";
			else GeneStateContent += stringtemp.str() + "]";
			stringtemp.clear();
		}
	}
	else
	{
		for (int x=1; x<101; x++)
		{
			git = find(GeneTypes->begin(), GeneTypes->end(), x);
			if (git != GeneTypes->end())	//Then we have this type in our genome.
			{
				index = distance(GeneTypes->begin(), git);
				git = GeneStates->begin();
				advance(git, index);
				stringstream stringtemp;
				stringtemp << (*git);
				if((size_t)count_types < GeneStates->size()-1) GeneStateContent += stringtemp.str() + ", ";
				else GeneStateContent += stringtemp.str() + "]";
				stringtemp.clear();
				count_types++;
			}
		}
	}
	return GeneStateContent;
}

string Genome::PrintGeneTypeContent()
{
	string GeneTypeContent = "[";
	for (int i=0; (size_t)i<GeneTypes->size(); i++)
	{
		stringstream stringtemp;
		stringtemp << GeneTypes->at(i);
		if((size_t)i<GeneTypes->size()-1) GeneTypeContent += stringtemp.str() + ", ";
		else GeneTypeContent += stringtemp.str() + "]";
		stringtemp.clear();
	}
	return GeneTypeContent;
}
