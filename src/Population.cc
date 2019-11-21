#include "Population.hh"

Population::Population() {
	int i,j;
	p_nr_proks_=0;

	for(i=0;i<NR;i++) for(j=0;j<NC;j++){
			PPSpace[i][j]=NULL;
	}
}

Population::~Population() {
	int i,j;
	for(i=0;i<NR;i++) for(j=0;j<NC;j++){
			if((PPSpace[i][j])!=NULL){
				delete (PPSpace[i][j]);
			}
	}
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
		PP_Copy=new Prokaryote();
		PP_Copy->ClonePPFromPP(PP);
		PPSpace[row][col] = PP_Copy;
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
		CP = new Prokaryote();
		CP->Replicate(PP);
		if(Parent_Genome != CP->G->PrintContent(NULL, false)) cout << "Child #" << n << ":\t" << CP->G->PrintContent(NULL, false) << endl;
		delete CP;
		CP = NULL;
	}
}

void Population::UpdatePopulation()	//This is the main next-state function.
{

	if(Time%TimeSaveGrid==0)	PrintSampleToFile();
	//if(Time%TimeSaveGrid==0)	PrintFieldToFile();
	if(Time%TimeTerminalOutput==0)	ShowGeneralProgress();

	for(int i=0; i<NR; i++) for(int j=0; j<NC; j++)		//Here we flag all individuals that are eligible for replication. Getting to the M-stage in this timestep only makes you eligible for replication in the next timestep.
	{
		if(PPSpace[i][j] != NULL)
		{
			if(PPSpace[i][j]->Stage == 4)	PPSpace[i][j]->ready_for_replication = true;
			else	PPSpace[i][j]->ready_for_replication = false;
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
				if (PPSpace[nrow][ncol]->ready_for_replication && chance < repl_rate)	//Only previously flagged individuals get to replicate (i.e. not the ones that acquired the M-stage only this timestep).
				{
					// cout << "Parent: " << PPSpace[nrow][ncol]->G->PrintGeneStateContent() << "\t" << PPSpace[nrow][ncol]->G->PrintContent(NULL, true) << endl;
					// cout << "Parent: " << PPSpace[nrow][ncol]->Stage << endl;
					PPSpace[i][j] = new Prokaryote();
					PPSpace[i][j]->Replicate(PPSpace[nrow][ncol]);
					// cout << "Parent: " << PPSpace[nrow][ncol]->Stage << endl;
					// cout << "Child: " << PPSpace[i][j]->Stage << endl;
					// cout << "Child: " << PPSpace[i][j]->G->PrintGeneStateContent() << "\t" << PPSpace[i][j]->G->PrintContent(NULL, true) << endl;
				}
			}

		}
		else	//Site is alive.
		{
			if(chance < death_rate)	//Death events
			{
				delete PPSpace[i][j];
				PPSpace[i][j]=NULL;
			}
			else	//Update InternalState of prokaryote
			{
				PPSpace[i][j]->G->UpdateGeneStates();
				PPSpace[i][j]->UpdateCellCycle();
			}
		}
	}

	//Do some diffusion here?

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
