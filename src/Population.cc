#include "Population.hh"

Population::Population() {
	int i,j;
	p_nr_proks_=0;
	p_id_count_=0;

	Fossils = new FossilRecord();
	for(i=0;i<NR;i++) for(j=0;j<NC;j++){
			PPSpace[i][j]=NULL;
	}
}

Population::~Population() {
	int i,j;
	for(i=0;i<NR;i++) for(j=0;j<NC;j++){
			if((PPSpace[i][j])!=NULL){
				delete (PPSpace[i][j]);
				PPSpace[i][j]=NULL;
			}
	}
	delete Fossils;
	Fossils=NULL;
}

void Population::InitialisePopulation()
{
	Prokaryote* PP;
	Prokaryote* PP_Copy;

	//First create one Prokaryote.
	PP = new Prokaryote();
	PP->InitialiseProkaryote();
	//Print its expression.
	cout << "\nInitial expression = " << PP->G->PrintGeneStateContent() << endl;
	//Print its genome.
	cout << "Initial genome = " << PP->G->PrintContent(NULL, true) << endl;	//Pass true as argument to print the nicely colored format for the terminal output.
	//Now fill the field with this prokaryote (I guess this is less intensive then creating new randomized prokaryotes for the whole grid).
	for(int row=0; row<NR; row++) for(int col=0; col<NC; col++){
		p_id_count_++;	//Make sure the first individual gets p_id_count_ of 1.
		PP_Copy=new Prokaryote();
		PP_Copy->ClonePPFromPP(PP, p_id_count_);
		PP_Copy->Ancestor = NULL;	//Null-pointer tells me the cell was initialised.
		PPSpace[row][col] = PP_Copy;
		Fossils->BuryFossil(PPSpace[row][col]);
		if(p_id_count_<=generation_sample)	OldGeneration[p_id_count_-1] = PPSpace[row][col];	//Put a subset of prokaryote pointers in the OldGeneration array.
	}

	delete PP;	//I cannot delete PP_Copy, because each is actually turned into one of grid spaces. I can however delete this single bit of memory.
	PP = NULL;
}

void Population::ReproduceMasterGenome()
{
	Prokaryote* PP;	//the parent
	Prokaryote* CP;	//the child
	//First create one Prokaryote.
	PP = new Prokaryote();
	PP->InitialiseProkaryote();
	//Print the genome to let the user be aware.
	string Parent_Genome = PP->G->PrintContent(NULL, false);
	cout << "\nMaster genome:\t" << Parent_Genome << endl;
	//Now make children one by one, through replication (incl. of course mutation) and print each to the terminal, before freeing the memory and reusing the pointer.
	//Maybe this is not very efficient, since the programme will be writing between each reproduction event; but it is definitely simple.
	for(int n=0; n<NrMutants; n++)
	{
		PP->Replicate();
		CP = new Prokaryote();
		CP->Mitosis(PP, 1);	//It does not matter which fossil_id each child gets, so I say 1.
		if(Parent_Genome != CP->G->PrintContent(NULL, false, false)) cout << "Child #" << n << ":\t" << CP->G->PrintContent(NULL, false, false) << endl;
		delete CP;
		CP = NULL;
	}
}

void Population::UpdatePopulation()	//This is the main next-state function.
{

	if(Time%TimeSaveGrid==0)	PrintSampleToFile();
	//if(Time%TimeSaveGrid==0)	PrintFieldToFile();
void Population::UpdatePopulation()	//This is the main next-state function.
{
	int reps = 0;
	if(Time%TimeSaveGrid==0)
	{
		if(NR*NC > 3000)	PrintSampleToFile();
		else	PrintFieldToFile();
	}
	if(Time%TimeTerminalOutput==0)	ShowGeneralProgress();
	if(Time%TimePruneFossils==0 && Time!=0)	PruneFossilRecord();
	if(Time%TimeOutputFossils==0 && Time!=0)	Fossils->ExhibitFossils();

	for(int i=0; i<NR; i++) for(int j=0; j<NC; j++)		//Here we flag all individuals that are eligible for replication. Getting to the M-stage in this timestep only makes you eligible for replication in the next timestep.
	{
		if(PPSpace[i][j] != NULL)
		{
			if(PPSpace[i][j]->Stage == 4)	PPSpace[i][j]->ready_for_division = true;
			else	PPSpace[i][j]->ready_for_division = false;
		}
	}

	int update_order[NR*NC];
	for(int u=0; u<NR*NC; u++) update_order[u]=u;
	random_shuffle(&update_order[0], &update_order[NR*NC]);		//Is also set by initial_seed through srand(initial_seed); see World.cc

	for(int u=0; u<NR*NC; u++)		//Go through the field: birth, death (and later possibly diffusion).
	{
		int i = update_order[u]/NR;	//Row index.
		int j = update_order[u]%NR;	//Column index.
		double chance = uniform();

		if (PPSpace[i][j]==NULL)	//Site is empty.
		{
			int random_neighbour = (int)(uniform()*9);
			int ni = random_neighbour/replication_neighbourhood;
			int nj = random_neighbour%replication_neighbourhood;

			//Wrap grid boundaries
			int nrow = i+ni-1;
			if(nrow < 0)	nrow = NR-1;
			else if(nrow >= NR)	nrow = 0;

			int ncol = j+nj-1;
			if(ncol < 0)	ncol = NC-1;
			else if(ncol >= NC)	ncol = 0;

			if (PPSpace[nrow][ncol]!=NULL)	//Random neighbour is alive.
			{
				//Replication events
				if (PPSpace[nrow][ncol]->ready_for_division && chance < (repl_rate<PPSpace[nrow][ncol]->fitness_deficit?0:1)*pow(repl_rate - PPSpace[nrow][ncol]->fitness_deficit,2))	//Only previously flagged individuals get to replicate (i.e. not the ones that acquired the M-stage only this timestep).
				{
					// cout << "Parent (2n): " << PPSpace[nrow][ncol]->G->PrintGeneStateContent() << "\t" << PPSpace[nrow][ncol]->G->PrintContent(NULL, false, false) << endl;
					PPSpace[i][j] = new Prokaryote();
					p_id_count_++;
					PPSpace[i][j]->Mitosis(PPSpace[nrow][ncol], p_id_count_);
					// cout << "Parent (n): " << PPSpace[nrow][ncol]->G->PrintGeneStateContent() << "\t" << PPSpace[nrow][ncol]->G->PrintContent(NULL, false, false) << endl;
					// cout << "Child (n): " << PPSpace[i][j]->G->PrintGeneStateContent() << "\t" << PPSpace[i][j]->G->PrintContent(NULL, false, false) << endl;
					if(PPSpace[i][j]->mutant)	Fossils->BuryFossil(PPSpace[i][j]);
					reps++;
				}
			}

		}
		else	//Site is alive.
		{
			if(chance < death_rate)	//Death events
			{
				if(!PPSpace[i][j]->mutant && !PPSpace[i][j]->saved_in_graveyard)
				{
					delete PPSpace[i][j];
				}
				else	PPSpace[i][j]->alive = false;
				PPSpace[i][j] = NULL;
			}

			else	//Update internal state of prokaryote
			{
				PPSpace[i][j]->G->UpdateGeneStates();
				PPSpace[i][j]->UpdateCellCycle();

				if (PPSpace[i][j]->Stage == 2)
				{
					PPSpace[i][j]->Replicate();
				}
			}
		}

	}
	//Do some diffusion here?
}

void Population::PruneFossilRecord()
{
	std::list<int> AllFossilIDs;
	/* How the record is pruned:
	*
	* MRCA is a pointer to an agent in the 'fossil record'. It points to the specific agent that started a new genotype (so only
	* mutants are stored in this ancestortrace.
	*
	* As long as the MRCA is not one of the first 50 generated (which have NULL as a common ancestor), keep on looking for
	* ancestors recursively. Add all IDs found like this to a big list (AllFossilIDs). Later, all individuals that are
	* not in the list, but are part of the fossil record, are deleted from the fossil record. Individuals that
	* are still alive are never deleted, since it is still unknown if they will be succesfull. (also see example asci)
	*

	BEFORE PRUNING:
	Ancestor list contains
	1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21


	PRUNING METHOD:
	Tracing back [living individuals] from t=x to t=0


	[ t=0 ]-->-->-->-->-->-->-->-->--[ t=x ]

	  2
	 /
	1		            14<<<<[ 15 ]
	 <<		          <<
	   3--4--5        10<<<<12<<<<<<<<[ 13 ]
	    <<	        <<
	       6<<<7<<<9
		\	<<
	        8         11<<<<16<<<<<<<<[ 17 ]
			   <<
			     18--19
			       <<
				 20<<<<<<<[ 21 ]

				 _________________
				|-- extinct branch|
				|<< traced branch |
				 -----------------

	AFTER PRUNING:
	Trace did not include:
	19,8,4,5,2
	Pruned ancestor list contains
	1,3,6,7,9,10,11,12,13,14,15,16,17,18,20,21

	 */
	for(int i=0; i<NR; i++)	for(int j=0; j<NC; j++)
	{
		if(PPSpace[i][j]!=NULL)
		{
			// Last common ancestor of living individual located
			Prokaryote* lastCA = PPSpace[i][j]->Ancestor;
			// Value 1 is always the root, which has a Null-parent.
			while(lastCA != NULL)
			{
				AllFossilIDs.push_back(lastCA->fossil_id);	// Added to list of agents we must keep
				lastCA = lastCA->Ancestor;	// 'Next level of taxonomy' (e.g. from 14 to 12 in example asci)
			}
		}
	}
	// Delete duplicates (e.g. Agent 9 in example asci will be located 4 times. Agent 11 two times, etc.)
	AllFossilIDs.sort();
	AllFossilIDs.unique();


	// Delete all in FossilList that are not in AllFossilIDs (unless they are still living):
	cout << "ID count: " << p_id_count_ << endl;
	cout << "Before pruning: " << (*Fossils).FossilList.size() << endl;
	iterpps ip = (*Fossils).FossilList.begin();
	while(ip != (*Fossils).FossilList.end())
	{
		// Search if stored agent was also found by tracing back:
		int fossilID = (*ip)->fossil_id;
		iter findit = std::find(AllFossilIDs.begin(),AllFossilIDs.end(),fossilID);
		// If not, delete the fossil unless it is still alive or is still saved in the graveyard. If a prokaryote dies, the graveyard-flag remains for one ShowGeneralProgress() cycle at most, so that the fossil can be deleted at the next pruning step. If ShowGeneralProgress() precedes PruneFossilRecord(), this is issue is even avoided, because flags are already removed off dead prokaryotes.
		if(findit==AllFossilIDs.end() && !(*ip)->alive && !(*ip)->saved_in_graveyard)
		{
			delete *ip;
			ip = Fossils->FossilList.erase(ip);
		}
		else
		{
			++ip;
		}
	}
	AllFossilIDs.clear();
	cout << "After pruning: " << (*Fossils).FossilList.size() << endl;
}

}

void Population::PrintFieldToFile()
{
	FILE* f;
	char OutputFile[800];
	sprintf(OutputFile, "%s/snapgrids/field%08d.txt", folder.c_str(), Time);
	f=fopen(OutputFile, "w");
	if (f == NULL){	printf("Failed to open file for writing the field.\n");	}

	for (int i=0; i<NR; i++) for(int j=0; j<NC; j++) {	//Don't print row and col numbers to save memory, these can be extracted by secondary scripts.
		if(PPSpace[i][j]==NULL){
		 	fprintf(f, "0\n");
		}
		else{		//Print internal state and genome of prokaryote to file.
			fprintf(f, "%s\t", PPSpace[i][j]->G->PrintGeneStateContent().c_str());
			fprintf(f, "%s\t", PPSpace[i][j]->G->PrintGeneTypeContent().c_str());
			fprintf(f, "%s\n", PPSpace[i][j]->G->PrintContent(NULL, false).c_str());
		}
	}
	fclose(f);
}

void Population::PrintSampleToFile()
{
	FILE* f;
	char OutputFile[800];
	sprintf(OutputFile, "%s/snapsamples/sample%08d.txt", folder.c_str(), Time);
	f=fopen(OutputFile, "w");
	if (f == NULL){	printf("Failed to open file for writing the sample.\n");	}

	int save_number = 3000;
	int count_saved = 0;
	for (int i=0; i<NR*2; i++) for(int j=0; j<NC; j++) {	//Don't print row and col numbers to save memory, these can be extracted by secondary scripts.
		if (uniform() < ((double)save_number/(double)(NR*NC)))
		{
			if(PPSpace[i%NR][j]==NULL){
			 	fprintf(f, "0\n");
			}
			else{		//Print internal state and genome of prokaryote to file.
				fprintf(f, "%s\t", PPSpace[i%NR][j]->G->PrintGeneStateContent().c_str());
				fprintf(f, "%s\t", PPSpace[i%NR][j]->G->PrintGeneTypeContent().c_str());
				fprintf(f, "%s\n", PPSpace[i%NR][j]->G->PrintContent(NULL, false).c_str());
			}
			count_saved ++;
			if(count_saved==save_number)
			{
				fclose(f);
				return;
			}
		}
	}

}

void Population::ShowGeneralProgress() {
	int alive=0;//, g1=0, s=0, g2=0, m=0, d=0;
	int stages[5] = {0, 0, 0, 0, 0};
	for (int i=0; i<NR; i++) for(int j=0; j<NC; j++) {
		if (PPSpace[i][j]!=NULL){
			alive++;
			stages[PPSpace[i][j]->Stage]++;
		}
	}
	printf("T %d\t() %d\t\tD %d\tG1 %d\tS %d\tG2 %d\tM %d\n", Time, alive, stages[0], stages[1], stages[2], stages[3], stages[4]);
	if (alive==0)
	{
		printf("And since there is no more life, we will stop the experiment here.\n\n");
		exit(1);
	}
}
