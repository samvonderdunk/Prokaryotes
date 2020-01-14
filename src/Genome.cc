#include "Genome.hh"

Genome::Genome() {
	BeadList=NULL;
	GeneStates=NULL;
	GeneTypes=NULL;
	gnr_genes=0;
	g_length=0;
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
		if(IsGene(*ii))
		{
			Gene* gene_template = dynamic_cast<Gene*>(*ii);
			Gene* new_gene = dynamic_cast<Gene*>(bead);
			new_gene->original_five = gene_template->original_five;		//It seems that this works: by making a gene-pointer to the bead, we change the bead, and transmit the original_five element.
		}
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

void Genome::CloneGenome(const Genome* G_template){
	BeadList=new list<Bead*>();		//The pointer ChromBBList (indeed of type std::list<ChromBB*>) is now initialised.
	CopyPartOfGenome(G_template->BeadList->begin(),G_template->BeadList->end());	//c is a pointer to a genome that we are going to clone. By dereferencing (->) we take from the actual Genome object (i.e. not its pointer), the ChromBBList object, which is itself a pointer to a list of pointers to beads (i.e. the main genome structure). Because ChromBBList is only a pointer to the list of bead-pointers, we dereference it again (->) to get the actual list. For the list we can then ask where the first element is (begin()) and when the list ends (end()).
	GeneStates=new vector<int>(G_template->GeneStates->size(),0);
	GeneTypes=new vector<int>(G_template->GeneStates->size(),0);
	for(int i=0; (size_t)i<G_template->GeneStates->size(); i++)
	{
		GeneStates->at(i) = G_template->GeneStates->at(i);
		GeneTypes->at(i) = G_template->GeneTypes->at(i);
	}

	gnr_genes = G_template->gnr_genes;
	g_length = G_template->g_length;
	mutated_bitstring = false;
	mutated_type = false;
	mutated_new_tfbs = false;
}

bool Genome::ReplicateGenome(){	//My own function that should double all the genes in the genome.
	list<Bead*>* Temp_BeadList;
	Temp_BeadList=new list<Bead*>();

	CopyPartOfGenomeToTemplate(BeadList->begin(), BeadList->end(), Temp_BeadList);

	//Append template to whole genome
	Bead *bead;
	iter ii = Temp_BeadList->begin();
	while(ii != Temp_BeadList->end()){
		bead=(*ii)->Clone();
		(*BeadList).push_back(bead);
		g_length++;
		ii++;
	}

	//Release memory
	iter i;
	i = Temp_BeadList->begin();
	while(i != Temp_BeadList->end()) {
		delete (*i);
		i++;
	}
	i = Temp_BeadList->erase(Temp_BeadList->begin(),Temp_BeadList->end());
	delete Temp_BeadList;
	Temp_BeadList=NULL;

	return true;
}

/*
###########################################################################
###########################################################################
							 |\ /\  |  | ----  /\  ---- o / \  |\ |
							|  V  \ l_J   |   /- \  |   | L_J  | \|
###########################################################################
###########################################################################
*/

void Genome::MutateGenome(){
	iter it;
	it=(*BeadList).begin();
	while(it!=(*BeadList).end()) {
		if(IsGene(*it)) {
			it=GeneMutate(it);
		}
		else if(IsTFBS(*it)) {
			it=TFBSMutate(it);
		}
		else {
			it++;
		}
	}
	//TFBSInnovation();	//We will have to define a similar thing for GeneInnovation later...
	if (g_length > 200)
	{
		printf("Warning: genome sizes reaching extravagant size (%d).\nExiting just to be safe...\n", g_length);
		exit(1);
	}

	if(mutated_bitstring || mutated_type) SetClaimVectors();	//In this case, a gene type may have been lost or gained, so the ClaimVectors need to be resized. Note that if a tfbs is duplicated its ClaimVector does not need to change and for tfbs deletions the ClaimVector is automatically deleted as well... (I think).
}

Genome::iter Genome::GeneMutate(iter ii) {
	Gene* gene;
	gene = dynamic_cast<Gene*>(*ii);
	bool potential_type_change = false;

	double uu = uniform();
	if(uu < gene_duplication_mu)
	{
		ii = GeneDuplication(ii);
		return ii;
	}

	else if(uu < gene_deletion_mu+gene_duplication_mu)
	{
		ii = GeneDeletion(ii);
		return ii;
	}

	else
	{
		if(uniform() < gene_threshold_mu)	//All mutational events are now independent, i.e. a gene can get multiple mutations or none, or one.
		{
			if (uniform()>0.8)	gene->threshold = (uniform()>0.5) ? gene->threshold+1 : gene->threshold-1;
			else	gene->threshold = (int)(uniform()*(2*WeightRange+1) - WeightRange);
		}

		if (uniform() < gene_activity_mu)
		{
			potential_type_change = true;
			if (uniform()>0.8){
				gene->activity = (uniform()>0.5) ? gene->activity+1 : gene->activity-1;
			}
			else{
				gene->activity = (uniform()>0.5) ? -1*uniform()*WeightRange : 1*uniform()*WeightRange;
			}
			potential_type_change = true;
		}

		for(int k=0; k<binding_length; k++)
		{
			if(uniform() < gene_binding_domain_mu)
			{
				if (gene->binding_domain[k] == false) gene->binding_domain[k] = true;
				else if (gene->binding_domain[k] == true) gene->binding_domain[k] = false;
			 	potential_type_change = true;
				mutated_bitstring = true;
			}
		}

		if (potential_type_change)
		{
			PotentialTypeChange(ii);// && !gene->original_five)	//Check that it has not become the same type as one of the other genes. I don't think it matters that it is one of the original five. Every gene can convert into an existing type; for duplicated genes, one can attain a new type.
		}
		ii++;
		return ii;
	}
}

Genome::iter Genome::TFBSMutate(iter ii)
{
	TFBS* tfbs;
	tfbs = dynamic_cast<TFBS*>(*ii);		//downcast the Bead-pointer to a TFBS-pointer with the name tfbs;

	double uu = uniform();
	if(uu < tfbs_duplication_mu)
	{
		ii = TFBSDuplication(ii);
		return ii;
	}
	else if(uu < tfbs_duplication_mu+tfbs_deletion_mu)
	{
		ii = TFBSDeletion(ii);
		return ii;
	}
	else
	{
		for (int k=0; k<binding_length; k++)
		{
			if(uniform() < tfbs_binding_site_mu)
			{
				if (tfbs->binding_site[k] == false) tfbs->binding_site[k] = true;
				else if (tfbs->binding_site[k] == true) tfbs->binding_site[k] = false;
				mutated_bitstring = true;
			}
		}
		if(uniform() < tfbs_activity_mu)
		{
			if (uniform()>0.8)	tfbs->activity = (uniform()>0.5) ? tfbs->activity+1 : tfbs->activity-1;
			else	tfbs->activity = (uniform()>0.5) ? -1*uniform()*WeightRange : 1*uniform()*WeightRange;
		}

		ii++;
		return ii;
	}
}

Genome::iter Genome::GeneDuplication(iter ii)	//Currently duplication will not directly increment the expression of the gene, but it should in the next update.
{
	int copy_length;
	iter insertsite, first, last;
	list<Bead*> BeadListTemp;	//Create a new temporary genome list.
	last = ii;
	last++;   //One further than the gene position (the one not to be part of the dupl).
	first=FindFirstTFBSInFrontOfGene(ii);	//first tfbs in front of gene
	copy_length=distance(first, last);		//function from std:: takes two iterators (but not one iterator and one reverse_iterator!).
	CopyPartOfGenomeToTemplate(first, last, &BeadListTemp); //Makes a 'chromosome' with only this gene (and its tfbs) on it
	insertsite=FindRandomGenePosition();			// find position of gene to insert in front of
	insertsite=FindFirstTFBSInFrontOfGene(insertsite);	// find fist tfbs in front of this gene
	BeadList->splice(insertsite, BeadListTemp);	// splice temporary list into chromosome

	ii=last;	// make sure ii points to one position further than just duplicated gene
	g_length+=copy_length;
	gnr_genes++;
	return ii;
}

Genome::iter Genome::GeneDeletion(iter ii)
{
	iter first, last, jj;
	last=ii;//gene position
	last++;//one further than the gene position
	first=FindFirstTFBSInFrontOfGene(ii);//first tfbs in front of gene
	g_length -= distance(first, last);

	int type_abundance = CountTypeAbundance((*ii)->type);
	if(type_abundance < 2) LoseGeneType((*ii)->type);	//If gene deletion means losing a gene type, shrink GeneStates and GeneTypes vectors.
	jj=first;
	while( jj != last )
	{
		delete *jj;		//Here you delete the element that the iterator is pointing to.
		jj++;
	}
	gnr_genes--;//you know one gene is removed
	ii=(*BeadList).erase(first, last);		//Here you remove the iterators, i.e. the pointers to the deleted elements from the iterator 'list'.
	return ii;
}

Genome::iter Genome::TFBSDuplication(iter ii)
{
	TFBS* tfbs;
	tfbs=dynamic_cast<TFBS *>(*ii);
	TFBS* tfbsnew = new TFBS(*tfbs);

	iter tt;
	iter upstream;
	int randpos;

	tt = (*BeadList).begin();
	upstream = (*BeadList).begin();
	randpos = (int)(uniform()*g_length);
	advance(tt,randpos);			// tt holds random spot in the genome e.g. |-------x--------------------|
	upstream = tt;
	if(tt != BeadList->begin()) upstream--;	//Making upstream point to the FRONT of the first bead-pointer? e.g. u|-------x--------------------------|

	tt = (*BeadList).insert(tt, tfbsnew);		//While I understand the insert function, I don't see why you would want tt to become a particular value at this point; probably it doesn't matter.
	g_length++;

	ii++;
	return ii;
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

Genome::iter Genome::FindFirstTFBSInFrontOfGene(iter ii) const
{
	reviter rii(ii);	//BEHIND in this case means an element to the left in the list (i.e if ii points to the 6th element, rii(ii) will make rii point to the 5th element). The function .base() as used below will make rii.base() point to the 6th element again.

	reviter jj = (*BeadList).rend();//search should be bounded
	while(rii != jj)//begin not yet reached
	{
		if(!IsGene(*rii))
		{
			rii++;
		}
		else
			jj = rii;
	}
	return jj.base();
}

Genome::iter Genome::FindRandomGenePosition() const
{
	std::list< iter > pos;
	std::list< iter >::iterator ipos;
	iter i;
	iter ii;
	int randpos;

	if(gnr_genes==0)	return (*BeadList).end();
	else
	{
		i=(*BeadList).begin();
		while(i != (*BeadList).end())
		{
			if(IsGene(*i))	pos.push_back(i);
			i++;
		}
		randpos=(int)(uniform()*gnr_genes);
		if(randpos > gnr_genes)
		{
			printf("Error: Random gene outside genome limits.\n");
			exit(1);
		}
		ii = (*boost::next(pos.begin(),randpos));	//Possibly try advance(ii, randpos) instead.

		return ii;
	}
}

void Genome::PotentialTypeChange(iter ii)
{
	Gene* gene = dynamic_cast<Gene*>(*ii);
	int type_abundance = CountTypeAbundance(gene->type);
	bool found_matching_type = false;
	iter jj = BeadList->begin();	//The other genes in the genome
	while(jj != BeadList->end())
	{
		if(IsGene(*jj) && jj!=ii)	//We don't convert genes to themselves.
		{
			bool genes_are_the_same = CheckSameGeneTypes(ii, jj);

			if(genes_are_the_same)
			{
				if(type_abundance < 2)	LoseGeneType(gene->type);
				gene->type = (*jj)->type;		//Convert to existing gene type.
				found_matching_type = true;
				return;		//We have found a match, converted the gene; time to try mutation of the next bead.
			}
		}
		jj++;
	}
	if(found_matching_type == false && type_abundance > 1)	//We haven't been able to convert it to an existing type, so let's define it as a new type || Or it remains the same type.
	{
		gene_iter git;
		for(int x=1; x<1000; x++)		//Find a not yet used number to use as the type.
		{
			git = find(GeneTypes->begin(), GeneTypes->end(), x);
			if(git == GeneTypes->end())	//X is not yet used as a gene type.
			{
				GeneStates->push_back(0);		//If gene is active it stays active after a mutation.
				gene->type = x;
				GeneTypes->push_back(x);
				mutated_type = true;		//Not actually, but ClaimVectors should definitely be updated.
				break;
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

void Genome::LoseGeneType(int type)	//Check whether a type is lost.
{
	gene_iter git = find(GeneTypes->begin(), GeneTypes->end(), type);
	int index = distance(GeneTypes->begin(), git);
	GeneTypes->erase(git);

	git = GeneStates->begin();
	advance(git, index);
	GeneStates->erase(git);
	mutated_type = true;		//This means that we well have to run SetClaimVectors() afterwards.
}

int Genome::CountTypeAbundance(int type)
{
	int count_type = 0;
	iter ii = BeadList->begin();
	while (ii != BeadList->end())
	{
		if(IsGene(*ii))
		{
			Gene* gene = dynamic_cast<Gene*>(*ii);
			if(gene->type==type)	count_type++;
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

	SetClaimVectors();
}


void Genome::ReadBeadsFromString(string genome)
{
	char* bead;
	int index;
	bool StageInit[5];
	gnr_genes = 0;
	g_length = 0;
	Gene* gene;
	TFBS* tfbs;
	int type, threshold, activity, q;
	char* buffer;
	bool bitstring[binding_length];

	if (genestate_init == "G1")	for(int g=1; g<6; g++)	StageInit[g] = StageTargets[0][g-1];
	else if (genestate_init == "S") for(int g=1; g<6; g++)	StageInit[g] = StageTargets[1][g-1];
	else if (genestate_init == "G2") for(int g=1; g<6; g++)	StageInit[g] = StageTargets[2][g-1];
	else if (genestate_init == "M") for(int g=1; g<6; g++)	StageInit[g] = StageTargets[3][g-1];
	else	printf("Expression not set to a cell-cycle stage.\n");

	bead = strtok((char*)genome.c_str(),".");
	while (bead != NULL)
	{
		if(bead[1] == 'G')	//Bead is a gene
		{
			buffer = new char();
			int success = sscanf(bead, "(G%d:%d:%d:%s)", &type, &threshold, &activity, buffer);
			if(success != 4) cerr << "Could not find sufficient information for this gene. Genome file potentially corrupt. \n" << endl;
			q = 0;
			while(buffer[q] != ')')		//The extra bracket stored in buffer actually pays off here because I don't know how else I would know that we reached the end of buffer.
			{
				bitstring[q] = (buffer[q]=='1');	//Easiest way I could think of to convert a character to a boolean; (bool)buffer[q] always returns 1 (whether '1' or '0')!
				q++;
			}
			gene = new Gene(type, threshold, activity, bitstring, 0);	//We initialise with zero expression, but these should be updated with the first round of UpdateGeneStates().
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
						GeneStates->push_back((int)(uniform()*2));	//For the genes other than the five main TFs we randomly set the expression.
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
			delete buffer;	//I think this prevents a memory leak. Otherwise buffer simply frees up new memory the next time it says "buffer = new char()". Because I did not use "new []" to create the variable, I should not use "delete []" to free up the memory.
		}
		else	//Bead is a tfbs
		{
			buffer = new char();
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
			delete buffer;
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
	g_length = 0;

	Gene* gene;
	TFBS* tfbs;
	int tfbs_type=0, activity, threshold;		//"tfbs_type=0": Binding is defined by bitstring so type does not matter anymore.
	bool original_five;
	bool binding_section[binding_length];

	//For now, I don't mind that genomes are ordered by gene type.
	for (int gene_type=1; gene_type<=init_nr_gene_types; gene_type++)
	{

		//Define its tfbs's.
		for (int j=0; j<init_nr_tfbs_per_gene; j++)
		{
			activity = (uniform()>0.5) ? -1 : 1;
			for (int k=0; k<binding_length; k++)	binding_section[k] = (uniform()>0.5) ? true : false;
			tfbs = new TFBS(tfbs_type, activity, binding_section);
			(*BeadList).push_back(tfbs);
			g_length++;
		}

		//Define the gene itself. Here the type matters, and this is the only case where we first want to say that it is a particular type and then make sure that it has a unique bitstring. This bitstring may mutate, in which case we should track it to know whether it is one of the five main types. Conversely, I could say that the five main types are always defined by five particular genes. If these genes the definition of the type also changes. Other genes will be of different types until they happen to randomly attain the same properties as one of the original five genes, in which case they are set to that type. All this is important for duplication and deletion events.
		threshold = (int)(uniform()*(2*WeightRange+1) - (int)WeightRange);	//Threshold between -3 and 3 (including these borders).
		activity = (uniform()>0.5) ? -1 : 1;
		for (int k=0; k<binding_length; k++)	binding_section[k] = (uniform()>0.5) ? true : false;
		original_five = (abs(gene_type) < 6) ? true: false;	//The first five types become the original_five.
		//Check that we have not initiated a gene of the same type
		iter jj = BeadList->begin();
		while(jj != BeadList->end())
		{
			if(IsGene(*jj))
			{
				Gene* check_gene = dynamic_cast<Gene*>(*jj);
				bool genes_are_the_same = true;		//Prove me wrong!
				if(activity != check_gene->activity) genes_are_the_same = false;
				for (int k=0; k<binding_length; k++)
				{
					if(binding_section[k] != check_gene->binding_domain[k]) genes_are_the_same = false;
				}

				if(genes_are_the_same)
				{
					printf("Had to remake G%d because it is the same as G%d!\n", gene_type, check_gene->type);
					for(int t=0; t<init_nr_tfbs_per_gene; t++)	(*BeadList).pop_back();
					g_length -= init_nr_gene_types;
					gene_type--;
					break;
				}
			}
			jj++;
		}

		//gene = new Gene(gene_type, threshold, activity, gene_binding_domain, expression, original_five);	//We should only get here if the new gene is not the same as one of the old genes.
		gene = new Gene(gene_type, threshold, activity, binding_section, original_five);
		(*BeadList).push_back(gene);
		gnr_genes++;
		g_length++;

	}
	//Randomly initialise gene expression [0,g] for each type (where g is the number of gene of that type).
	GeneStates = new vector<int>();
	GeneTypes = new vector<int>();
	for(int g=1; g<=init_nr_gene_types; g++)
	{
		GeneStates->push_back((int)(uniform()*2));
		GeneTypes->push_back(g);
	}

	//For each TFBS, put all bindings strengths with all genes in its claim vector.
	SetClaimVectors();
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
			//Find the right index for GeneTypes/GeneStates of the gene type that we have in our hands.
			gene_iter git = find(GeneTypes->begin(), GeneTypes->end(), gene->type);
			int index = distance(GeneTypes->begin(), git);
			NextStateExpression->at(index) += max(min(sumeffects+1,1),0);	//For the +1, see explanation of the gene threshold in Gene.hh
			//(*NextStateExpression).push_back(max(min(sumeffects,1),0));		//For now we limit expression of individual genes to 1 or 0, i.e. on or off.
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
	double binding_strength;
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
					Gene* gene = dynamic_cast<Gene*>(*ig);
					gene_iter git = find(GeneTypes->begin(), GeneTypes->end(), gene->type);
					int index = distance(GeneTypes->begin(), git);
					binding_strength = MatchBitStrings(*it, *ig);
					tfbs->ClaimVector->at(index) = binding_strength;
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
	double claim_sum = 0.0;
	iter i_gene;
	//Calculate the claim sum.	//Maybe faster to iterate through the elements.
	for(int g=0; (size_t)g<GeneStates->size(); g++)
	{
		claim_sum += tfbs->ClaimVector->at(g) * (double) GeneStates->at(g);	//Claim of gene types is proportional to the number of active genes of that type and to their binding strength to the tfbs.
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
			if ( die_roll <= tfbs->ClaimVector->at(g) * (double)(GeneStates->at(g)) )
			{
				//Increment i_gene until we find a gene of the correct type.
				while(!(IsGene(*i_gene) && (*i_gene)->type==GeneTypes->at(g)))
				{
					i_gene++;
				}
				return i_gene;
			}
			die_roll -= tfbs->ClaimVector->at(g) * (double)GeneStates->at(g);
		}
	}
	printf("Error: claim out of bounds.\n");
	exit(1);
}

double Genome::MatchBitStrings(Bead* b_tfbs, Bead* b_gene){
	TFBS* tfbs;
	Gene* gene;
	tfbs = dynamic_cast<TFBS*> (b_tfbs);
	gene = dynamic_cast<Gene*> (b_gene);

  int binding_strength = 0;
  for(int k=0; k<binding_length; k++)
	{
    if(tfbs->binding_site[k] != gene->binding_domain[k]) binding_strength++;		//We do complementary pairing for the beauty of it (but also expected in visualisations).
  }

  return pow((double)binding_strength/(double)binding_length, tfbs_selection_exponent);
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
	int index, type, gene_act, gene_copynum, bead_count = 0;
	TFBS* tfbs;
	TFBS::claim_iter cit;
	double Effect;
	Gene* gene;

	//Before we can interpret the ClaimVectors of all TFBS's we need to know if there are genes present in higher copy number.
	TypeCopyNumber.resize(GeneTypes->size(),0);
	GeneActivity.resize(GeneTypes->size(),0);

	it = BeadList->begin();
	while(it != BeadList->end())
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
	while(it != BeadList->end())
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
				Effect = (*cit) * gene_copynum * gene_act * tfbs->activity;	//(*cit) is the raw claim, i.e. the bitstring match between TFBS and gene.

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

int Genome::FindIndexOfType(int type)
{
	gene_iter git = find(GeneTypes->begin(), GeneTypes->end(), type);
	if(git == GeneTypes->end())	return -1;	//This is a signal during MoveGenomeToChild() that the type does not yet exist in the GeneTypes vector.
	else	return distance(GeneTypes->begin(), git);
}

string Genome::PrintContent(list<Bead*> * chromosome, bool terminal)	// Printing function for terminal output.
{
	string GenomeContent;
	if(chromosome == NULL) chromosome = this->BeadList;
	iter i = chromosome->begin();
	Gene *gene;
	TFBS *tfbs;

	string gene_color_prefix = "";
	string gene_color_suffix = "";
	string tfbs_color_prefix = "";
	string tfbs_color_suffix = "";

	if(terminal==true){
		gene_color_prefix = "\033[94m";
		gene_color_suffix = "\033[0m";
		tfbs_color_prefix = "\033[92m";
		tfbs_color_suffix = "\033[0m";
	}

	while(i!=chromosome->end()) {
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

		GenomeContent += ")";
		i++;
	}
	return GenomeContent;
}

string Genome::PrintGeneStateContent()
{
	string GeneStateContent = "[";
	for (int i=0; (size_t)i<GeneStates->size(); i++)
	{
		stringstream stringtemp;
		stringtemp << GeneStates->at(i);
		if((size_t)i<GeneStates->size()-1) GeneStateContent += stringtemp.str() + ", ";
		else GeneStateContent += stringtemp.str() + "]";
		stringtemp.clear();
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
