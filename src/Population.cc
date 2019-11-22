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

void Population::ContinuePopulationFromBackup()
{
	ReadBackupFile();
	ReadAncestorFile();	//Currently, the fossil_ids are missing from the backup-file so it is impossible to link the fossils to live prokaryotes. But in the new version this will be possible.
	PruneFossilRecord();
}

void Population::ReadBackupFile()
{
	ifstream infile(backup_reboot.c_str());
	string line, data;
	char* data_element;
	string::iterator sit;
	Genome::iter it;
	int reading = 1;	//0 - nothing, 1 - GeneStates, 2 - GeneTypes, 3 - Genome, 4 - Prokaryote properties, 5 - Gene expression.
	int read_integer = 0, index, begin_data, end_data, counter, success, stage, pfork, panti_ori;
	bool is_mutant, is_mutant_child;
	unsigned long long prok_id;
	double deficit;
	Prokaryote* PP;
	Gene* gene;

	if (!infile.is_open())
	{
		printf("Backup-file could not be opened.\n");
		exit(1);
	}

	printf("Reading backup from file: %s\n", backup_reboot.c_str());
	int count_lines = 0;
	while(getline(infile,line))
	{
		if (count_lines/NC >= NR || count_lines%NR >= NC)
		{
			printf("Backup file was larger than the field; aborting just to be safe.\n");
			exit(1);
		}
		if(line == "0")	PPSpace[count_lines/NC][count_lines%NR] = NULL;
		else
		{
			//Start new individual.
			PP = new Prokaryote();
			PP->EmptyProkaryote();
			PP->G->BeadList = new list<Bead*>();
			PP->G->GeneTypes = new vector<int>();
			PP->G->GeneStates = new vector<int>();

			//Read BeadList. Do this first, because it first sets GeneStates and GeneTypes (randomised or based on other input files). We can reset these below.
			begin_data = line.find_first_of("(");
			end_data = line.find_last_of(")");
			data = line.substr(begin_data, end_data-begin_data+1);
			PP->G->ReadBeadsFromString(data);

			//Read GeneStates.
			index = line.find("]");
			data = line.substr(1, index-1);
			counter = 0;
			sit = data.begin();
			while(sit != data.end())
			{
				if(*sit == ',')
				{
					PP->G->GeneStates->push_back(read_integer);	//Save value just read.
					counter++;
					read_integer = 0;
				}
				else if(*sit != ' ')
				{
					read_integer *= 10;
					read_integer += (int)*sit - 48;
				}
				sit++;
			}
			PP->G->GeneStates->push_back(read_integer);
			read_integer = 0;

			//Read GeneTypes.
			begin_data = line.find("[", index);
			end_data = line.find("]", index+1);
			data = line.substr(begin_data+1, end_data-begin_data-1);
			counter = 0;
			sit = data.begin();
			while(sit != data.end())
			{
				if(*sit == ',')
				{
					PP->G->GeneTypes->at(counter) = read_integer;	//Save value just read.
					counter++;
					read_integer = 0;
				}
				else if(*sit != ' ')
				{
					read_integer *= 10;
					read_integer += (int)*sit - 48;
				}
				sit++;
			}
			PP->G->GeneTypes->at(counter) = read_integer;
			read_integer = 0;

			//Read prokaryote data.
			begin_data = line.find_last_of("[");
			end_data = line.find_last_of("]");
			data = line.substr(begin_data+1, end_data-begin_data-1);
			data_element = strtok((char*)data.c_str(),"\t");
			while(data_element != NULL)
			{
				success = sscanf(data_element, "%d %lf %d %d %llu %d %d", &stage, &deficit, &pfork, &panti_ori, &prok_id, &is_mutant, &is_mutant_child);
				if(success != 7)
				{
					cerr << "Could not find sufficient information for this prokaryote. Backup file potentially corrupt.\n" << endl;
					exit(1);
				}
				data_element = strtok(NULL, "\t");
				PP->Stage = stage;
				PP->fitness_deficit = deficit;
				PP->G->pos_fork = pfork;
				PP->G->pos_anti_ori = panti_ori;
				PP->fossil_id = prok_id;
				PP->mutant = is_mutant;
				PP->mutant_child = is_mutant_child;
				if (prok_id > p_id_count_) p_id_count_ = prok_id;
			}

			// Read expression of individual genes. This will work in the new version, where the nr of genes actually corresponds to the length of the gene expression data.
			begin_data = line.find("{");
			end_data = line.find("}");
			data = line.substr(begin_data+1, end_data-begin_data-2);
			sit = data.begin();
			it = PP->G->BeadList->begin();
			while (sit != data.end())
			{
				if(*sit != ' ')
				{
					while(!PP->G->IsGene(*it))	it++;	//Go through beads until you hit the next gene.
					gene = dynamic_cast<Gene*>(*it);
					gene->expression = (int)*sit - 48;
				}
				sit++;
			}

			if(PP->Stage == 2)	//Not sure if this is the right condition.
			{
				PP->G->MutationList = new vector<bool>(PP->G->pos_anti_ori);	//If you initiate MutationList during the programme (i.e. first time you get to ReplicateGenomeStep()), you will actually just make the MutationList the same length as g_length. But now we are reading in prokaryotes that have already replicated some beads, so that there g_length is longer and does not match the MutationList data in the backup file.

				begin_data = line.find_last_of("{");
				end_data = line.find_last_of("}");
				data = line.substr(begin_data+1, end_data-begin_data-3);	//The -3 is very strange (see how it was -2 above..) but seems to work now.
				sit = data.begin();
				counter = 0;
				read_integer = 0;

				while (sit != data.end())
				{
					if (*sit == ' ')
					{
						if(counter == 0)	PP->G->deletion_length = read_integer;
						else	PP->G->MutationList->at(counter-1) = (read_integer==1) ? true:false;
						read_integer = 0;
						counter++;
					}
					else	//We are looking at a number supposedly.
					{
						read_integer *= 10;
						read_integer += (int)*sit - 48;
					}
					sit++;
				}

				PP->G->MutationList->at(counter-1) = (read_integer==1) ? true:false;
				read_integer = 0;

			}

			PP->G->SetClaimVectors();

			PPSpace[count_lines/NC][count_lines%NR] = PP;
			if (PP->mutant)	Fossils->BuryFossil(PP);
		}
		if(count_lines<generation_sample)
		{
			if(PPSpace[count_lines/NC][count_lines%NR] != NULL)	PPSpace[count_lines/NC][count_lines%NR]->saved_in_graveyard = true;
			OldGeneration[count_lines] = PPSpace[count_lines/NC][count_lines%NR];	//Put a subset of prokaryote pointers in the OldGeneration array.
		}
		count_lines++;
	}
	if (count_lines!= NR*NC)
	{
		printf("Backup file was too small for the field; aborting just to be safe.\n");
		exit(1);
	}
}

void Population::ReadAncestorFile()
{
	ifstream infile(anctrace_reboot.c_str());
	string line, data;
	int begin_data, end_data, TimeOA;
	unsigned long long ID, AncID;
	iterpps ip, ip2;
	Prokaryote* PP;

	if (!infile.is_open())
	{
		printf("Ancestor-file could not be opened.\n");
		exit(1);
	}

	printf("Reading ancestors from file: %s\n", anctrace_reboot.c_str());
	int count_lines = 0;
	int count_alive = 0;
	int count_fossils = 0;
	while(getline(infile,line))
	{
		count_lines++;
		end_data = line.find("\t");
		data = line.substr(0,end_data);
		stringstream(data) >> ID;

		begin_data = end_data;
		end_data = line.find("\t",end_data+1);
		data = line.substr(begin_data, end_data-begin_data);
		stringstream(data) >> AncID;

		begin_data = end_data;
		end_data = line.find("\t",end_data+1);
		data = line.substr(begin_data, end_data-begin_data);
		stringstream(data) >> TimeOA;

		begin_data = end_data;
		end_data = line.size();
		data = line.substr(begin_data+1, end_data-begin_data);

		ip = Fossils->FossilList.begin();
		while (ip != Fossils->FossilList.end())
		{
			if ((*ip)->fossil_id == ID)	//Then we have found a live prokaryote in our ancestor file.
			{
				count_alive++;
				(*ip)->time_of_appearance = TimeOA;
				if (AncID == 0)	(*ip)->Ancestor = NULL;
				else
				{
					ip2 = Fossils->FossilList.begin();
					while(ip2 != Fossils->FossilList.end())
					{
						if ((*ip2)->fossil_id == AncID)
						{
							(*ip)->Ancestor = *ip2;
							break;
						}
						ip2++;
					}
					if(ip2 == Fossils->FossilList.end())
					{
						printf("Error: ancestor not found...exiting.\n");
						exit(1);
					}
				}
				break;
			}
			ip++;
		}

		if (ip == Fossils->FossilList.end())	//We did not break out of the loop, so we have apparently not encountered this ID among the current list of fossils.
		{
			count_fossils++;
			PP = new Prokaryote();
			PP->EmptyProkaryote();
			PP->fossil_id = ID;
			PP->alive = false;
			PP->mutant = true;
			PP->time_of_appearance = TimeOA;
			if(AncID == 0)	PP->Ancestor = NULL;
			else	//We have to find the rightful parent of this creature. It should be in the list because we read in the fossil record starting with the oldest fossils; everything afterwards should have a parent present in the record.
			{
				ip2 = Fossils->FossilList.begin();
				while (ip2 != Fossils->FossilList.end())
				{
					if ((*ip2)->fossil_id == AncID)
					{
						PP->Ancestor = *ip2;
						break;
					}
					ip2++;
				}
				if(ip2 == Fossils->FossilList.end())
				{
					printf("Error: ancestor not found...exiting.\n");
					exit(1);
				}
			}

			//Set up its ghost genome.	It only has beads.
			PP->G->BeadList = new list<Bead*>();
			PP->G->GeneTypes = new vector<int>();
			PP->G->GeneStates = new vector<int>();
			PP->G->ReadBeadsFromString(data);
			Fossils->BuryFossil(PP);
		}

	}
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
	if(Time%TimeSaveBackup==0 && Time!=0)	OutputBackup();

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
	sprintf(OutputFile, "%s/snapsamples/field%08d.txt", folder.c_str(), Time);
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
	for (int i=0; i<NR*2; i++) for(int j=0; j<NC; j++) {
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

void Population::OutputBackup()
{
	Genome::iter it;
	Genome::mut_iter mit;
	Gene* gene;
	FILE* f;
	char OutputFile[800];
	sprintf(OutputFile, "%s/backups/backup%08d.txt", folder.c_str(), Time);
	f=fopen(OutputFile, "w");
	if (f == NULL)	printf("Failed to open file for writing the backup.\n");

	for (int i=0; i<NR; i++) for(int j=0; j<NC; j++) {
		if(PPSpace[i][j]==NULL){
		 	fprintf(f, "0\n");
		}
		else	//Print internal state and genome of prokaryote to file.
		{
			fprintf(f, "%s\t", PPSpace[i][j]->G->PrintGeneStateContent().c_str());
			fprintf(f, "%s\t", PPSpace[i][j]->G->PrintGeneTypeContent().c_str());
			fprintf(f, "%s\t", PPSpace[i][j]->G->PrintContent(NULL, false, false).c_str());
			fprintf(f, "[%d %f %d %d %llu %d %d]\t", PPSpace[i][j]->Stage, PPSpace[i][j]->fitness_deficit, PPSpace[i][j]->G->pos_fork, PPSpace[i][j]->G->pos_anti_ori, PPSpace[i][j]->fossil_id, PPSpace[i][j]->mutant, PPSpace[i][j]->mutant_child);
			fprintf(f, "{");
			it = PPSpace[i][j]->G->BeadList->begin();
			while (it != PPSpace[i][j]->G->BeadList->end())
			{
				if(PPSpace[i][j]->G->IsGene(*it))
				{
					gene = dynamic_cast<Gene*>(*it);
					fprintf(f, "%d ", gene->expression);
				}
				it++;
			}
			if(PPSpace[i][j]->G->MutationList == NULL)	fprintf(f, "\b}\n");
			else
			{
				fprintf(f, "\b}\t{%d ", PPSpace[i][j]->G->deletion_length);
				mit = PPSpace[i][j]->G->MutationList->begin();
				while (mit != PPSpace[i][j]->G->MutationList->end())
				{
					fprintf(f, "%d ", ((*mit)==true?1:0));
					mit++;
				}
				fprintf(f, "\b}\n");
			}
		}
	}

	fclose(f);
	return;
}

void Population::ShowGeneralProgress()
{
	int alive=0, live_comparisons=0, present_alives=0;//, g1=0, s=0, g2=0, m=0, d=0;
	int stages[5] = {0, 0, 0, 0, 0};
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
