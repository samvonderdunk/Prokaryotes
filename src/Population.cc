#include "Population.hh"

Population::Population()
{
	int i,j,n;
	p_nr_proks_=0;
	p_id_count_=0;

	Fossils = new FossilRecord();
	for(i=0;i<NR;i++) for(j=0;j<NC;j++)
	{
			PPSpace[i][j]=NULL;
	}
	for(n=0;n<generation_sample;n++)
	{
		OldGeneration[n]=NULL;
	}
	Environment = init_env;
}

Population::~Population()
{
	int i,j,n;
	iterpps ips;

	for(i=0;i<NR;i++) for(j=0;j<NC;j++)
	{
		if((PPSpace[i][j])!=NULL)
		{
			if(PPSpace[i][j]->saved_in_graveyard)	//With diffusion I have to find the location in the OldGeneration array, bc it is not per se the same (i.e. cannot be determined from the same i, j).
			{
				for (n=0;n<generation_sample;n++)
				{
					if (OldGeneration[n]==PPSpace[i][j])
					{
						OldGeneration[n]=NULL;
						break;
					}
				}
			}
			if(PPSpace[i][j]->mutant)	Fossils->EraseFossil(PPSpace[i][j]->fossil_id);
			delete (PPSpace[i][j]);
			PPSpace[i][j]=NULL;
		}
	}

	for(n=0;n<generation_sample;n++)
	{
		if(OldGeneration[n]!=NULL)
		{
			if(OldGeneration[n]->mutant)	Fossils->EraseFossil(OldGeneration[n]->fossil_id);
			delete (OldGeneration[n]);
			OldGeneration[n]=NULL;
		}
	}

	delete Fossils;	//This will delete the rest of FossilList in the Fossils class internally.
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
	cout << "\nInitial expression = " << PP->G->PrintGeneStateContent(false) << endl;
	//Print its genome.
	cout << "Initial genome = " << PP->G->PrintContent(NULL, true, false) << endl;	//Pass true as argument to print the nicely colored format for the terminal output.
	//Now fill the field with this prokaryote (I guess this is less intensive then creating new randomized prokaryotes for the whole grid).
	for(int row=0; row<NR; row++) for(int col=0; col<NC; col++){
		// if(row<20 && col<20)	//Initialise square
		if (uniform() < 0.1)	//Initialise lower density
		{
			p_id_count_++;	//Make sure the first individual gets p_id_count_ of 1.
			PP_Copy=new Prokaryote();
			PP_Copy->ClonePPFromPP(PP, p_id_count_);
			PP_Copy->Ancestor = NULL;	//Null-pointer tells me the cell was initialised.
			PPSpace[row][col] = PP_Copy;
			Fossils->BuryFossil(PPSpace[row][col]);
			if (row*NC+col <= generation_sample)
			{
				OldGeneration[row*NC+col] = PPSpace[row][col];	//Put a subset of prokaryote pointers in the OldGeneration array.
				PPSpace[row][col]->saved_in_graveyard = true;
			}
		}
	}

	delete PP;	//I cannot delete PP_Copy, because each is actually turned into one of grid spaces. I can however delete this single bit of memory.
	PP = NULL;

	if (environmental_noise)	NoiseEnvironment();
	cout << "Initial environment = " << Environment << endl;
}


//Soon, priviliges will have to be read from the backup-file as well.
void Population::ContinuePopulationFromBackup()
{
	ReadBackupFile();
	if(anctrace_reboot != "")	ReadAncestorFile();	//Currently, the fossil_ids are missing from the backup-file so it is impossible to link the fossils to live prokaryotes. But in the new version this will be possible.
	else	//Otherwise we will reset the fossil_ids to prevent these getting out of hand.
	{
		p_id_count_ = 0;
		for (int i=0; i<NR; i++)	for(int j=0; j<NC; j++)
		{
			if (PPSpace[i][j]!=NULL)
			{
				p_id_count_++;
				PPSpace[i][j]->fossil_id = p_id_count_;	//Other things such as time_of_appearance and Ancestor are set to zero by the EmptyProkaryote function.
			}
		}
	}

	OutputBackup();
	PruneFossilRecord();

	if (environmental_noise)	NoiseEnvironment();
}

void Population::ReadBackupFile()
{
	ifstream infile(backup_reboot.c_str());
	string line, data;
	char* data_element;
	string::iterator sit;
	Genome::iter it;
	int i, read_integer = 0, index, begin_data, end_data, counter, success, stage, pfork, panti_ori, temp_is_mutant, temp_priv, nr=NR, nc=NC, init_seed, read_header = 0;	//For now, set nr and nc for backup-file to NR and NC parameters (i.e. if backup-file does not contain header; for old backups). May be removed later.
	bool is_mutant, priv;
	unsigned long long prok_id, sdraws;
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

		//Read "Header" of backup file.
		if (line=="### Header ###")
		{
			read_header = 1;
			continue;
		}
		else if (line=="#### Main ####")
		{
			read_header = -1;
			continue;
		}
		else if (read_header==1)
		{
			data_element = (char*)line.c_str();
			success = sscanf(data_element, "NR:%d\tNC:%d", &nr, &nc);
			if(success != 2)
			{
				cerr << "Could not read NR and NC from backup-file.\n" << endl;
				exit(1);
			}
			else	cout << "nr=" << nr << ", nc=" << nc << endl;
			read_header = 2;
			continue;
		}
		else if (read_header==2)
		{
			data_element = (char*)line.c_str();
			success = sscanf(data_element, "Initial seed:%d\tSeed draws:%llu", &init_seed, &sdraws);
			if (success != 2)
			{
				cerr << "Could not read seed status from backup-file.\n" << endl;
				exit(1);
			}
			else
			{
				cout << "Simulating " << sdraws << " random draws, initial seed=" << init_seed << endl;
				for(i=0; (unsigned long long)i<sdraws; i++){
					uniform();
				}
			}
			read_header = 3;	//Not doing anything with this yet.
			continue;
		}

		//Start reading "Main" data from backup file.
		//Make sure that field dimensions of backup and current simulation are related to one another.
		if (count_lines/nc >= NR || count_lines%nc >= NC)
		{
			count_lines++;
			continue;		//Sites that lie on [x,y] in the backup-sim, but where either x or y is not within the limits of the current sim; we simply lose these individuals.
		}
		if (count_lines%10000 == 0 && count_lines!=0)	printf("%d\n", count_lines);	//Print some progress.
		if(line == "0")	PPSpace[count_lines/nc][count_lines%nc] = NULL;	//Empty sites.
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
				success = sscanf(data_element, "%d %lf %d %d %llu %d %d", &stage, &deficit, &pfork, &panti_ori, &prok_id, &temp_is_mutant, &temp_priv);	//WARNING: in newer backup-files, the AncID has been inserted before temp_is_mutant and temp_priv, so BEWARE!!!
				if(success != 7)
				{
					cerr << "Could not find sufficient information for this prokaryote. Backup file potentially corrupt.\n" << endl;
					exit(1);
				}
				is_mutant = temp_is_mutant;	//Integers to bools (maybe not completely necessary..?)
				priv = temp_priv;

				data_element = strtok(NULL, "\t");
				PP->Stage = stage;
				PP->fitness_deficit = deficit;
				PP->G->pos_fork = pfork;
				PP->G->pos_anti_ori = panti_ori;
				PP->fossil_id = prok_id;
				PP->mutant = is_mutant;
				PP->priviliges = priv;
				if (PP->Stage>=2)	PP->time_replicated = 1;	//NOTE: we're not doing anything with time_replicated currently, so it's fine to set everything to 1 (allowing mitosis).
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
					while(!PP->G->IsGene(*it))
					{
						it++;	//Go through beads until you hit the next gene.
					}

					gene = dynamic_cast<Gene*>(*it);
					gene->expression = (int)*sit - 48;
				}
				else	it++;
				sit++;
			}

			PP->G->SetClaimVectors();

			PPSpace[count_lines/nc][count_lines%nc] = PP;
			if (PP->mutant)	Fossils->BuryFossil(PP);
		}
		if(count_lines<generation_sample)
		{
			if(PPSpace[count_lines/nc][count_lines%nc] != NULL)	PPSpace[count_lines/nc][count_lines%nc]->saved_in_graveyard = true;
			OldGeneration[count_lines] = PPSpace[count_lines/nc][count_lines%nc];	//Put a subset of prokaryote pointers in the OldGeneration array.
		}
		count_lines++;
	}
	if (count_lines!= nr*nc)
	{
		cout << count_lines << " " << nr << " " << nc << endl;
		printf("Length of backup file conflicts with header (nr*nc).\n");
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
	int count_alive = 0;
	int count_fossils = 0;
	int count_lines = 0;
	while(getline(infile,line))
	{
		if (count_lines%10000 == 0 && count_lines!=0)	printf("%d\n", count_lines);
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
			if ((*ip)->fossil_id == ID)	//Then we have found a live prokaryote in our ancestor file, because it will have to be added to the FossilRecord from the backup file. For these guys we only have to find its ancestor in the FossilList (all other data has been read from the backup-file).
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
						printf("Error: ancestor %llu not found...exiting.\n", ID);
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
					printf("Error: ancestor %llu not found...exiting.\n", ID);
					exit(1);
				}
			}

			//Set up its ghost genome.	It only has beads.
			PP->G->BeadList = new list<Bead*>();
			PP->G->GeneTypes = new vector<int>();
			PP->G->GeneStates = new vector<int>();
			PP->G->ReadBeadsFromString(data);
			PP->G->pos_anti_ori = PP->G->g_length;	//Otherwise the function PrintContent(NULL, false, true) (i.e. printing only the parental genome to a file) while think the parental genome is non-existent (pos_anti_or = 0).
			Fossils->BuryFossil(PP);
		}
		count_lines++;
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
	string Parent_Genome = PP->G->PrintContent(NULL, false, false);
	cout << "\nMaster genome:\t" << Parent_Genome << endl;
	//Now make children one by one, through replication (incl. of course mutation) and print each to the terminal, before freeing the memory and reusing the pointer.
	//Maybe this is not very efficient, since the programme will be writing between each reproduction event; but it is definitely simple.
	for(int n=0; n<NrMutants; n++)
	{
		while(PP->G->pos_fork < PP->G->pos_anti_ori)	PP->Replicate(100);	//One replication step should be enough for most genomes (i.e. GL <= 100).
		CP = new Prokaryote();
		CP->Mitosis(PP, 1);	//It does not matter which fossil_id each child gets, so I say 1.
		if(Parent_Genome != CP->G->PrintContent(NULL, false, false)) cout << "Child #" << n << ":\t" << CP->G->PrintContent(NULL, false, false) << endl;
		delete CP;
		CP = NULL;
	}
	delete PP;
	PP = NULL;
}

void Population::ScanMutationalPath()
{
	Prokaryote* IP;
	Prokaryote* PP;
	Prokaryote* CP;
	bool PhenotypicChange;	//Definition of "phenotype" can be one of two: 1) qualitative network configuration OR 2) expression pattern through time.
	int MutAttempts;	//If we have to try to many mutations to find one that maintains the phenotype we choose to stop.

	IP = new Prokaryote();
	IP->EmptyProkaryote();
	IP->InitialiseProkaryote();
	PP = new Prokaryote();
	PP->EmptyProkaryote();
	PP->G->CloneGenome(IP->G);
	cout << "\nInitial genome:\t" << IP->G->PrintContent(NULL, false, false) << endl;

	for(Time=0; Time<SimTime+1; Time++)	//We will do as many steps over the mutational path as given by SimTime.
	{
		CP = new Prokaryote();
		CP->EmptyProkaryote();
		CP->G->CloneGenome(PP->G);	//We start with the "child" being the same as the parent
		MutAttempts = 0;
		while( (CP->G->PrintContent(NULL, false, false) == PP->G->PrintContent(NULL, false, false)) || PhenotypicChange)	//We only select a new individual if it does not have the same genome as its parent (the current prokaryote) but is also not too different phenotypically.
		{
			delete CP;
			CP = NULL;
			while(PP->G->pos_fork < PP->G->pos_anti_ori) PP->Replicate(100);
			CP = new Prokaryote();
			CP->EmptyProkaryote();
			CP->Mitosis(PP, 1);
			// PhenotypicChange = QualitativeNetworkChange(IP, CP);	//One can also choose for the more dynamic option of comparing the child with its parent, but then the definition of the network might change over time.
			PhenotypicChange = CompareExpressionProgression(IP, CP, 4);	//Option 2: update both genomes' expression for x timesteps, starting in the expression provided through the command-line. If the second genome never deviated from the first in terms of the 5 main TFs, there has been no qualitative phenotypic change. You can compare over a larger number of timesteps if you want both to do many exactly synchronized cycles. This option requires a lot more trying, so runs slower.
			MutAttempts++;
			if(MutAttempts >= 1000000)
			{
				cout << "Stopped scanning mutational path -- path too narrow (" << MutAttempts << " attempts)" << endl;
				exit(1);
			}
		}
		cout << Time+1 << " steps:\t" << CP->G->PrintContent(NULL, false, false) << endl;
		//CP becomes the new "parent" for the next step, so we lose the data held by PP, make it the same as CP, and then erase CP, so that the next timestep we can make it a new child.
		delete PP;
		PP = NULL;
		PP = new Prokaryote();
		PP->EmptyProkaryote();
		PP->G->CloneGenome(CP->G);
		delete CP;
		CP = NULL;
	}
}

bool Population::CompareExpressionProgression(Prokaryote* PP1, Prokaryote* PP2, int time_window)
{
	int qualitative_change = false;
	Genome::gene_iter git;

	for (int wtime=0; wtime<time_window; wtime++)
	{
		PP1->G->UpdateGeneStates();
		PP2->G->UpdateGeneStates();

		//Check that the 2nd genome has the 5 main types. We assume the 1st genome does.
		for (int gt=1; gt<=5; gt++)
		{
			git = find(PP2->G->GeneTypes->begin(), PP2->G->GeneTypes->end(), gt);
			if (git == PP2->G->GeneTypes->end())
			{
				qualitative_change = true;
				break;
			}
		}

		//Only compare the first 5 genes' expression. For now it matters how many copies are expressed bc we stupidly compare these strings.
		if (PP1->G->PrintGeneStateContent(true).substr(1, 13) != PP2->G->PrintGeneStateContent(true).substr(1, 13))
		{
			qualitative_change = true;
			break;
		}
	}

	return qualitative_change;
}

void Population::FollowSingleIndividual()
{
	Prokaryote* PP;

	PP = new Prokaryote();
	PP->InitialiseProkaryote();

	for(Time=0; Time<SimTime+1; Time++)
	{
		if(environmental_noise)	NoiseEnvironment();
		//All we want is to know the expression pattern at each time step.
		cout << "T " << Time << "\tE " << Environment << "\tStage: " << PP->Stage << "\tPriviliges: " << PP->priviliges << "\tG_len: " << PP->G->g_length << "\tExpr: " << PP->G->PrintGeneStateContent(true) << endl;

		PP->G->UpdateGeneStates();
		if (PP->Stage <= 4)	PP->UpdateCellCycle();

		if (PP->Stage == 2 && PP->priviliges == true)
		{
			PP->Replicate(repl_step_size-Environment);
			PP->time_replicated++;
		}
		else if (PP->Stage > 4)	//Marked for immediate or future death.
		{
			cout << "T " << Time << "\tE " << Environment << "\tStage: -1\tPriviliges: " << PP->priviliges << "\tG_len: " << PP->G->g_length << "\tExpr: " << PP->G->PrintGeneStateContent(true) << endl;
			PP->Abortion();	//Although the cell is killed we can just start again by doing an abortion.
		}
		else if (PP->Stage == 4 && PP->priviliges == true && uniform() < repl_rate-PP->fitness_deficit && !(DivisionProtocol==2 && !uniform()<0.1))	//Attempting division. Since there are no other living beings, I don't have to consider the different DivisionProtocols or pick a random neighbour. Actually if we normally wait for an empty site, let's introduce a 10% chance to actually divide every timestep in M.
		{
			cout << "T " << Time << "\tE " << Environment << "\tStage: 4\tPriviliges: " << PP->priviliges << "\tG_len: " << PP->G->g_length << "\tExpr: " << PP->G->PrintGeneStateContent(true) << endl;
			PP->Abortion();	//The cell succesfully divides, but we don't need to create an actual child to continue looking at cell-cycle dynamics.
		}
	}

	delete PP;
	PP = NULL;
}

void Population::ExploreAttractorLandscape()
{
	int max_nr_initial_states, nr_initial_states, state_id, count_b, index;
	Genome::gene_iter git;
	Genome::iter it;
	Gene* gene;
	Prokaryote* PP;	//the parent
	Prokaryote* CP;	//the child
	vector<int> InitialState;
	vector< vector<int> > NextStates;
	vector< vector<int> >::iterator ns;
	vector<int> TimesVisited;
	bool novel_state;


	//First create one Prokaryote.
	PP = new Prokaryote();
	PP->InitialiseProkaryote();
	//Print the genome to let the user be aware.
	string Parent_Genome = PP->G->PrintContent(NULL, false, false);
	cout << "\nGenome:\t" << Parent_Genome << "\n" << endl;
	cout << "<<<\nSource state\tTarget state\tProbability" << endl;

	//Now run through all possible initial states (or sample the number of initial states given by the user).
	max_nr_initial_states = pow(2, PP->G->gnr_genes);
	nr_initial_states = (NrInitialStates > max_nr_initial_states) ? max_nr_initial_states : NrInitialStates;

	for (int s=0; s<nr_initial_states; s++)
	{
		//NextStates and TimesVisited are reset.
		NextStates.clear();
		TimesVisited.clear();

		//Set the GeneStates vector to zero, so that we only have to go through it when we need to turn a gene on (they are off by default).
		git = PP->G->GeneStates->begin();
		while(git != PP->G->GeneStates->end())
		{
			(*git) = 0;
			git++;
		}

		state_id = (NrInitialStates > max_nr_initial_states) ? s : (int)(uniform()*max_nr_initial_states);	//Sample states randomly (currently with redraws) or iterate through all possible states.

		//Translate state_id to a bit-pattern used to set the GeneStates.
		for(int b=PP->G->gnr_genes-1; b>=0; b--)
		{

			if (state_id >= pow(2, b))	//This gene is on.
			{
				//Find the b'th gene in the genome. We increment the expression of that type.
				count_b = -1;
				it = PP->G->BeadList->begin();
				while(it != PP->G->BeadList->end())
				{
					if(PP->G->IsGene(*it))
					{
						count_b++;
						if (count_b == b)	//We have hit the b'th gene.
						{
							git = find(PP->G->GeneTypes->begin(), PP->G->GeneTypes->end(), (*it)->type);
							index = distance(PP->G->GeneTypes->begin(), git);
							gene = dynamic_cast<Gene*>(*it);
							gene->expression = 1;
							PP->G->GeneStates->at(index) += gene->expression;
							break;	//Break out of the while-loop; we have set the expression of one gene according to the initial state that we want to test here.
						}
					}
					it++;
				}
				state_id -= pow(2, b);
			}

		}
		//Now we have set GeneStates based on the bit-pattern.
		InitialState = *PP->G->GeneStates;

		//Simulate this cell as many types as the grid is large.
		for(int a=0; a<NR*NC; a++)
		{
			CP=new Prokaryote();
			CP->ClonePPFromPP(PP, 0);
			CP->G->UpdateGeneStates();	//Update the expression.

			//Check if we have already seen this expression as an outcome.
			novel_state = true;	//We'll set it to false if we find out that we have already seen it.
			ns = NextStates.begin();
			while(ns != NextStates.end())
			{
				if (*(CP->G->GeneStates) == *ns)
				{
					index = distance(NextStates.begin(), ns);
					TimesVisited.at(index)++;
					novel_state = false;
					break;
				}
				ns++;
			}
			if (novel_state)	//Either it is the first state in general we have visited or it is a novel state.
			{
				NextStates.push_back(*CP->G->GeneStates);
				TimesVisited.push_back(1);
			}

			delete CP;
			CP = NULL;
		}

		//Now we have simulated this initial state NR*NC times; print the output.
		for (int p=0; (size_t)p<TimesVisited.size(); p++)
		{
			git = InitialState.begin();
			while (git != InitialState.end())
			{
				cout << *git;
				git++;
			}
			cout << "\t";
			git = (NextStates.at(p)).begin();
			while (git != (NextStates.at(p)).end())
			{
				cout << *git;
				git++;
			}
			cout << "\t" << (double)TimesVisited.at(p)/(double)(NR*NC) << endl;
		}

	}

	cout << ">>>\n" << endl;
	delete PP;
	PP = NULL;
}

/* ######################################################################## */
/*				UPDATE	UPDATE	UPDATE	UPDATE	UPDATE	UPDATE	UPDATE						*/
/* ######################################################################## */

void Population::UpdatePopulation()	//This is the main next-state function.
{
	int update_order[NR*NC];
	int u, i, j, e_block, s_block;
	double chance, resource, diffusion_steps;

	if(Time==TimeZero)	ResetProgressCounters();	//Initialise some of my output stats for the first time.

	if(Time%TimeSaveGrid==0)	PrintFieldToFile();
	if(Time%TimeTerminalOutput==0)	ShowGeneralProgress();
	if(Time%TimePruneFossils==0 && Time!=0)	PruneFossilRecord();
	if(Time%TimeOutputFossils==0 && Time!=0)	Fossils->ExhibitFossils();
	if(Time%TimeSaveBackup==0 && Time!=0)	OutputBackup();

	if(environmental_noise) NoiseEnvironment();	//Potential change of environment.

	ResetProgressCounters();

	for(u=0; u<NR*NC; u++) update_order[u]=u;
	random_shuffle(&update_order[0], &update_order[NR*NC]);		//Is also set by initial_seed through srand(initial_seed); see World.cc

	for(u=0; u<NR*NC; u++)		//Go through the field: birth, death (and later possibly diffusion).
	{
		i = update_order[u]/NC;	//Row index.
		j = update_order[u]%NC;	//Column index.
		e_block = j/(NC/env_blocks);
		s_block = j/(NC/stats_in_blocks);
		chance = uniform();

		Environment = gradient[e_block];

		if (PPSpace[i][j] != NULL)	//Alive site.	Programme is written such that nothing happens at empty sites.
		{
			if (chance < death_rate)
			{
				DeathOfProkaryote(i, j);
				continue;	//A death event ends the dynamics for this site at this point in time.
			}

			//Pick random neighbour (needed in some protocols to determine if cell is allowed to replicate).
			coords neigh = PickNeighbour(i, j);

			if (IsReadyToDivide(i, j, neigh.first, neigh.second) == true)	//Will the cell attempt division?
			{
				if (PPSpace[i][j]->Stage == 5)
				{
					nr_death_cycles[s_block]++;		//You cycled yourself to death (i.e. failed cycle).
					DeathOfProkaryote(i, j);	//Prokaryote was marked for death upon division (incomplete cycle).
					continue;
				}

				//If you get here, all is good, division is actually going to happen!
				if(PPSpace[i][j]->nr_offspring == 0)
				{
					cum_time_alive[s_block] += Time - PPSpace[i][j]->time_of_appearance;
					nr_first_births[s_block]++;
				}
				cum_fit_def[s_block] += PPSpace[i][j]->fitness_deficit;

				if (PPSpace[neigh.first][neigh.second] != NULL)	DeathOfProkaryote(neigh.first, neigh.second);	//Depending on the DivisionProtocol, we're committed to overgrowing neighbours.

				PPSpace[neigh.first][neigh.second] = new Prokaryote();
				p_id_count_++;
				PPSpace[neigh.first][neigh.second]->Mitosis(PPSpace[i][j], p_id_count_);
				if(PPSpace[neigh.first][neigh.second]->mutant)	Fossils->BuryFossil(PPSpace[neigh.first][neigh.second]);
				nr_birth_events[s_block]++;
			}

			PPSpace[i][j]->G->UpdateGeneStates();
			if (PPSpace[i][j]->Stage <= 4)	PPSpace[i][j]->UpdateCellCycle();	//Stage 5 is marked for death upon division (no second chances).

			if (PPSpace[i][j]->Stage == 2 && PPSpace[i][j]->priviliges == true)
			{
				resource = CollectResource(i, j, Environment);
				PPSpace[i][j]->Replicate(resource);
			}

			if (PPSpace[i][j]->Stage == 6)
			{
				nr_death_cycles[s_block]++;		//You cycled yourself to death (i.e. failed cycle).
				DeathOfProkaryote(i, j);	//Cell was marked for immediate death (i.e. no waiting for attempted division).
			}

		}
	}

	//Margolus diffusion.
	diffusion_steps = diffusion_rate;
	while (diffusion_steps > 0.)
	{
		if (diffusion_steps >= 1.)
		{
			MargolusDiffusion();
			diffusion_steps -= 1.;
		}
		else
		{
			if (uniform() < diffusion_steps)	MargolusDiffusion();
			diffusion_steps -= 1.;	//Make sure that it is not slightly above zero (?)
		}
	}

}

void Population::MargolusDiffusion()	//Based on Brem's function from Evolvabear_V9.0
{
	int kernel, i, j;
	Prokaryote* PPtemp;

	//Rotation of 2x2 squares divided up according to two different kernels (shift by one), and both run the entire field.
	//kernel=0 does not need to deal with wrapped boundaries if NR and NC are even, but kernel=1 does.
	for (kernel=0; kernel<2; kernel++)
	{
		for(i=kernel; i<NR-1+kernel; i+=2)	for(j=kernel; j<NC-1+kernel; j+=2)
		{
			if (uniform()<0.5)	//Counter-clockwise.
			{
				PPtemp = PPSpace[(i+1)%NR][j%NC];	//The %NR and %NC make sure that the grid is wrapped (e.g. the boundary row index "100" will return 100%100 = 0 if NR=100).
				PPSpace[(i+1)%NR][j%NC] = PPSpace[i%NR][j%NC];
				PPSpace[i%NR][j%NC] = PPSpace[i%NR][(j+1)%NC];
				PPSpace[i%NR][(j+1)%NC] = PPSpace[(i+1)%NR][(j+1)%NC];
				PPSpace[(i+1)%NR][(j+1)%NC] = PPtemp;
				PPtemp = NULL;
			}
			else								//Clockwise.
			{
				PPtemp = PPSpace[i%NR][(j+1)%NC];
				PPSpace[i%NR][(j+1)%NC] = PPSpace[i%NR][j%NC];
				PPSpace[i%NR][j%NC] = PPSpace[(i+1)%NR][j%NC];
				PPSpace[(i+1)%NR][j%NC] = PPSpace[(i+1)%NR][(j+1)%NC];
				PPSpace[(i+1)%NR][(j+1)%NC] = PPtemp;
				PPtemp = NULL;
			}
		}
	}
}


bool Population::IsReadyToDivide(int i, int j, int nrow, int ncol)	//Figure out if cell meets division criteria for current division protocol (see Header.hh).
{
		//Division is attempted in Stage 4 or higher (where 5 marks for death), and with probability determined by repl_rate corrected for fitness_deficit.
		//Note that the below statement is negative, so it uses "OR" constructions.
	if (PPSpace[i][j]->Stage < 4 || PPSpace[i][j]->priviliges == false || uniform() > (repl_rate - PPSpace[i][j]->fitness_deficit))
	{
		return false;
	}

	else	//Options are sorted from most to least aggressive.
	{
		if (DivisionProtocol == 0)	//Overgrow neighbour.
		{
			return true;
		}
		else if (DivisionProtocol == 1)	//Compete with neighbour.
		{
			if (PPSpace[nrow][ncol] == NULL) return true;
			else if (PPSpace[i][j]->time_stationary > PPSpace[nrow][ncol]->time_stationary)	return true;	//Competition rule.
			else
			{
				PPSpace[i][j]->Stage = 6;	//Failing in competition is lethal.
				return false;
			}
		}
		else if (DivisionProtocol == 2)	//Wait for empty site.
		{
			if (PPSpace[nrow][ncol] == NULL)	return true;
			else	return false;
		}
	}
}

void Population::DeathOfProkaryote(int i, int j)
{
	if(!PPSpace[i][j]->mutant && !PPSpace[i][j]->saved_in_graveyard)	delete PPSpace[i][j];
	else	PPSpace[i][j]->alive = false;
	PPSpace[i][j] = NULL;
}

void Population::NoiseEnvironment()
{
	if(uniform() < environmental_change_rate)	//Change environment.
	{
		Environment = (double)((int)(uniform()*(2*environmental_variation+1) - environmental_variation));
	}
}

double Population::CollectResource(int i, int j, double Environment)
{
	double resource, repl_noise;

	resource = repl_step_size - Environment;

	if (repl_step_noise)
	{
		repl_noise = ( (int)(uniform()*2) == 1 ) ? resource:-resource;
		resource = ( (int)(uniform()*2) == 1 ) ? resource:resource+repl_noise;
	}

	if (resource_dependent_replication)
	{
		resource *= 8 - NeighbourhoodDensity(i, j);
		resource /= 8;
	}

	return resource;
}

Population::coords Population::PickNeighbour(int i, int j)	//Pick random neighbour.
{
	int nrow = i, ncol = j, ni, nj, random_neighbour;

	while (nrow == i && ncol == j)	//Try again if you pick yourself.
	{
		random_neighbour = (int)(uniform()*9);
		ni = random_neighbour/replication_neighbourhood;
		nj = random_neighbour%replication_neighbourhood;

		//Wrap grid boundaries
		nrow = i+ni-1;
		if(nrow < 0)	nrow += NR;
		else if(nrow >= NR)	nrow -= NR;

		ncol = j+nj-1;
		if(ncol < 0)	ncol += NC;
		else if(ncol >= NC)	ncol -= NC;
	}
	return std::make_pair(nrow, ncol);
}

int Population::NeighbourhoodDensity(int i, int j)
{
	int ii, jj, nrow, ncol, density=0;

	for (ii=i-1; ii<=i+1; ii++) for (jj=j-1; jj<=j+1; jj++)
	{
		if (ii == i && jj == j)	continue;

		if (ii < 0)	nrow = ii + NR;	//-1 becomes -1+100=99, i.e. the last index of a row with 100 sites (0-99).
		else if (ii >= NR)	nrow = ii - NR;
		else	nrow = ii;
		if (jj < 0)	ncol = jj + NC;
		else if (jj >= NC)	ncol = jj - NC;
		else	ncol = jj;

		if (PPSpace[nrow][ncol] != NULL)	density++;
	}

	return density;
}

void Population::ResetProgressCounters()
{
	for(int i=0; i<stats_in_blocks; i++)
	{
		nr_birth_events[i] = 0;
		nr_death_cycles[i] = 0;
		cum_fit_def[i] = 0.0;
		nr_first_births[i] = 0;
		cum_time_alive[i] = 0;	//Store the total time that prokaryotes undergoing mitosis are alive (used to extract the average length of their life cycle).
	}

}

/* ######################################################################## */
/*				FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS	FOSSILS						*/
/* ######################################################################## */

void Population::PruneFossilRecord()
{
	std::list<unsigned long long> AllFossilIDs;
	int fossil_record_size;

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
	fossil_record_size = (*Fossils).FossilList.size();
	iterpps ip = (*Fossils).FossilList.begin();
	while(ip != (*Fossils).FossilList.end())
	{
		// Search if stored agent was also found by tracing back:
		unsigned long long fossilID = (*ip)->fossil_id;
		iterull findit = std::find(AllFossilIDs.begin(),AllFossilIDs.end(),fossilID);
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
	cout << "ID count (" << p_id_count_ << ")\tFossil record (pruned from " << fossil_record_size << " to " << (*Fossils).FossilList.size() << ")" << endl;
}

double Population::MatrixDistance(Prokaryote* PP1, Prokaryote* PP2)
{
	double** N1;
	double** N2;
	N1 = new double* [5];
	N2 = new double* [5];
	for (int row=0; row<5; row++)
	{
		N1[row] = new double[8];
		N2[row] = new double[8];
	}

	for(int row=0; row<5; row++)	for(int col=0; col<8; col++)
	{
		N1[row][col] = .0;
		N2[row][col] = .0;
	}
	PP1->G->GenomeToNetwork(N1);
	PP2->G->GenomeToNetwork(N2);

	double Distance = .0;
	for(int row=0; row<5; row++)	for(int col=0; col<8; col++)
	{
		Distance += abs(N1[row][col] - N2[row][col]);
	}

	for (int row=0; row<5; row++)
	{
		delete [] N1[row];
		N1[row] = NULL;
		delete [] N2[row];
		N2[row] = NULL;
	}
	delete [] N1;
	N1 = NULL;
	delete [] N2;
	N2 = NULL;

	return Distance;
}

bool Population::QualitativeNetworkChange(Prokaryote* PP1, Prokaryote* PP2)	//PP1 is used as a qualitative (predefined) network. If it has a positive interaction between two genes, PP2 should have an interaction that is at least +0.5 (-0.5 for negative, or between -0.5 and 0.5 for no interaction).
{
	bool qualitative_change = false;
	double** N1;
	double** N2;
	N1 = new double* [5];
	N2 = new double* [5];
	for (int row=0; row<5; row++)
	{
		N1[row] = new double[8];
		N2[row] = new double[8];
	}

	for(int row=0; row<5; row++)	for(int col=0; col<8; col++)
	{
		N1[row][col] = .0;
		N2[row][col] = .0;
	}
	PP1->G->GenomeToNetwork(N1);
	PP2->G->GenomeToNetwork(N2);

	for(int row=0; row<5; row++)	for(int col=0; col<8; col++)
	{
		if (N1[row][col] < -0.1 && N2[row][col] > -0.5)	qualitative_change = true;
		else if ( (N1[row][col] >= -0.1 && N1[row][col] <= 0.1) && (N2[row][col] <= -0.5 || N2[row][col] >= 0.5) ) qualitative_change = true;
		else if (N1[row][col] > 0.1 && N2[row][col] < 0.5) qualitative_change = true;
	}

	for (int row=0; row<5; row++)
	{
		delete [] N1[row];
		N1[row] = NULL;
		delete [] N2[row];
		N2[row] = NULL;
	}
	delete [] N1;
	N1 = NULL;
	delete [] N2;
	N2 = NULL;

	return qualitative_change;
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
		else{		//Print relevant variables to make snapshots.
			fprintf(f, "%d\t%d\t%d\t%d\t%d\n", PPSpace[i][j]->Stage, PPSpace[i][j]->G->g_length, PPSpace[i][j]->G->gnr_genes, PPSpace[i][j]->G->pos_anti_ori, PPSpace[i][j]->maturing_time);
		}
	}
	fclose(f);
}

void Population::PrintSampleToFile()	//Use this if the field is very large and you want to calculate field averages.
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
				fprintf(f, "%s\t", PPSpace[i%NR][j]->G->PrintGeneStateContent(false).c_str());
				fprintf(f, "%s\t", PPSpace[i%NR][j]->G->PrintGeneTypeContent().c_str());
				fprintf(f, "%s\n", PPSpace[i%NR][j]->G->PrintContent(NULL, false, true).c_str());
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
	Gene* gene;
	FILE* f;
	char OutputFile[800];
	unsigned long long AncestorID;

	sprintf(OutputFile, "%s/backups/backup%08d.txt", folder.c_str(), Time);
	f=fopen(OutputFile, "w");
	if (f == NULL)	printf("Failed to open file for writing the backup.\n");

	//Print NR and NC to file.
	fprintf(f, "### Header ###\nNR:%d\tNC:%d\nInitial seed:%d\tSeed draws:%llu\n#### Main ####\n", NR, NC, initial_seed, seed_draws);

	for (int i=0; i<NR; i++) for(int j=0; j<NC; j++) {
		if(PPSpace[i][j]==NULL){
		 	fprintf(f, "0\n");
		}
		else	//Print internal state and genome of prokaryote to file.
		{
			fprintf(f, "%s\t", PPSpace[i][j]->G->PrintGeneStateContent(false).c_str());
			fprintf(f, "%s\t", PPSpace[i][j]->G->PrintGeneTypeContent().c_str());
			fprintf(f, "%s\t", PPSpace[i][j]->G->PrintContent(NULL, false, false).c_str());

			if (PPSpace[i][j]->Ancestor==NULL)	AncestorID = 0;
			else	AncestorID = PPSpace[i][j]->Ancestor->fossil_id;
			fprintf(f, "[%d %f %d %d %llu %llu %d %d]\t", PPSpace[i][j]->Stage, PPSpace[i][j]->fitness_deficit, PPSpace[i][j]->G->pos_fork, PPSpace[i][j]->G->pos_anti_ori, PPSpace[i][j]->fossil_id, AncestorID, PPSpace[i][j]->mutant, PPSpace[i][j]->priviliges);
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
			fprintf(f, "\b}\n");
		}
	}

	fclose(f);
	return;
}

void Population::ShowGeneralProgress()
{
	int u, i, j, ii, jj, block, total_alive=0, live_comparisons=0, present_alives=0;
	int alive[stats_in_blocks];
	int stages[7][stats_in_blocks];
	string stage_names[7] = {"D", "G1", "S", "G2", "M", "D1", "D2"};
	double pop_distance = .0, pop_msd = .0;
	bool already_saved_in_graveyard;

	//Reset counters for total nr of proks and for each stage.
	for (i=0; i<stats_in_blocks; i++)
	{
		alive[i] = 0;
		for (j=0; j<7; j++)	stages[j][i] = 0;
	}

	for(i=0; i<NR; i++) for(j=0; j<NC; j++)
	{
		if (Time!=0 && (i*NC + j)<generation_sample)
		{
			if (OldGeneration[i*NC+j] != NULL)
			{
				if (PPSpace[i][j] != NULL)
				{
					//Compare current individual of grid with the one that was here at the previous measure frame.
					pop_distance += MatrixDistance(OldGeneration[i*NC+j], PPSpace[i][j]);
					live_comparisons++;
				}
				//OldGeneration[i*NR+j] is an actual prokaryote; before we overwrite we have to see whether it is time to completely delete this guy from all records (when it is not a mutant and not alive, thus not in the fossilrecord) or whether we just stop saving it in the graveyard (allowing the fossilrecord to decide whether it is still interesting for the geneology).
				if(!OldGeneration[i*NC+j]->mutant && !OldGeneration[i*NC+j]->alive)
				{
					delete OldGeneration[i*NC+j];	//If this fellow was dead, remove the grave if it is not interesting for the fossil record.
					OldGeneration[i*NC+j]=NULL;
				}
				else
				{
					//If it has moved somewhere else it might already be included in the next generation sample (bc it is alive at an "upstream" location), so then we should not set "saved_in_graveyard" to "false".
					already_saved_in_graveyard = false;
					for (u=0; u<i*NC+j; u++)
					{
						if (OldGeneration[i*NC+j] == PPSpace[u/NC][u%NC])
						{
							already_saved_in_graveyard = true;
							break;
						}
					}
					if (!already_saved_in_graveyard)	OldGeneration[i*NC+j]->saved_in_graveyard = false;	//We have to free them from this constrain; otherwise they can never be thrown out of the fossil record.
				}
			}
			OldGeneration[i*NC+j] = PPSpace[i][j];	//Update OldGeneration for the next ShowGeneralProgress. We also do this if one of these pointers was NULL.
			if(PPSpace[i][j] != NULL)	PPSpace[i][j]->saved_in_graveyard = true;	//Even if it dies, it will be kept around at least until the next generation has been compared to it (or longer if it is interesting for the fossil record).

			if ((i*NC + j) < pow(generation_sample, 0.5))	//In a 100x100 grid, the entire first row (100 individuals) will be compared amongst themselves.
			{
				for (ii=0; ii<NR; ii++)	for(jj=0; jj<NC; jj++)
				{
					if ((ii*NC + jj) < pow(generation_sample, 0.5))
					{
						if(PPSpace[i][j] != NULL && PPSpace[ii][jj] != NULL)
						{
							pop_msd += MatrixDistance(PPSpace[i][j], PPSpace[ii][jj]);
							present_alives++;
						}
					}
				}
			}
		}

		if (PPSpace[i][j]!=NULL)
		{
			block = j/(NC/stats_in_blocks);
			total_alive++;
			alive[block]++;
			stages[PPSpace[i][j]->Stage][block]++;
		}
	}

	//Printing. cout is used to achieve printing during running of the programme (vs. printf).
	cout << "T " << Time;		//Standard output, independent of number of blocks (Environment will be random in case there are multiple environments in the field, i.e. when there is a gradient).
	PrintStatPerBlock("E", gradient, NULL);	//Passing a double array. If stats_in_blocks is lower than the actual number of blocks in the gradient, only the first x environments from gradient will be printed.
	PrintStatPerBlock("()", alive, NULL);
	for (i=0; i<7; i++)	PrintStatPerBlock(stage_names[i], stages[i], NULL);
	PrintStatPerBlock("Births", nr_birth_events, NULL);
	PrintStatPerBlock("Death-cycs", nr_death_cycles, NULL);
	PrintStatPerBlock("FitD", cum_fit_def, nr_birth_events);	//Passing a double array.
	PrintStatPerBlock("CycleLen", cum_time_alive, nr_first_births);
	cout << "\tDist " << pop_distance/live_comparisons << "\tMSD " << pop_msd/present_alives << endl;

	if (total_alive==0)
	{
		cout << "And since there is no more life, we will stop the experiment here.\n" << endl;
		exit(1);
	}
}

void Population::PrintStatPerBlock(string stat_name, int* stat_array, int* scaler_array)
{
	int i;
	cout << "\t" << stat_name << " ";
	for (i=0; i<stats_in_blocks; i++)
	{
		if(scaler_array==NULL)	cout << *(stat_array+i);
		else	cout << (double)*(stat_array+i) / *(scaler_array+i);

		if (i<stats_in_blocks-1)	cout << ",";
	}
}

void Population::PrintStatPerBlock(string stat_name, const double* stat_array, int* scaler_array)
{
	int i;
	cout << "\t" << stat_name << " ";
	for (i=0; i<stats_in_blocks; i++)
	{
		if(scaler_array==NULL)	cout << *(stat_array+i);
		else	cout << (double)*(stat_array+i) / *(scaler_array+i);

		if (i<stats_in_blocks-1)	cout << ",";
	}
}
