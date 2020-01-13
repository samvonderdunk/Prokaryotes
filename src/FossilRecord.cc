#include "FossilRecord.hh"

FossilRecord::FossilRecord()
{
}

FossilRecord::~FossilRecord()
{
	iterpps ips;
	ips = FossilList.begin();
	while (ips != FossilList.end())
	{
		if (*ips != NULL)
		{
			delete (*ips);
			ips = FossilList.erase(ips);
		}
		else	ips++;
	}
}

void FossilRecord::EraseFossil(unsigned long long fossilID)
{
	iterpps ips;
	ips = FossilList.begin();
	while (ips != FossilList.end())
	{
		if ((*ips)->fossil_id == fossilID)
		{
			ips = FossilList.erase(ips);
			return;
		}
		ips++;
	}
}

void FossilRecord::BuryFossil(Prokaryote *P)
{
	FossilList.push_back(P);
}

void FossilRecord::ExhibitFossils()
{
	FILE* f;
	char OutputFile[800];
	sprintf(OutputFile, "%s/ancestors/anctrace%08d.txt", folder.c_str(), Time);
	f=fopen(OutputFile, "w");

	if (f == NULL)	printf("Failed to open file for writing the ancestor trace.\n");
	// fprintf(f, "#id\t#anc_id\t#time_oa\t#genome\n");	//Don't print the header to save space, but at least you know now!

	iterpps ip = FossilList.begin();
	while(ip != FossilList.end())
	{
		if ((*ip)->Ancestor == NULL)
		{
			fprintf(f, "%llu\t%d\t%d\t%s\n", (*ip)->fossil_id, 0, (*ip)->time_of_appearance, (*ip)->G->PrintContent(NULL, false, true).c_str());
		}
		else
		{
			fprintf(f, "%llu\t%llu\t%d\t%s\n", (*ip)->fossil_id, ((*ip)->Ancestor)->fossil_id, (*ip)->time_of_appearance, (*ip)->G->PrintContent(NULL, false, true).c_str());
		}
		ip++;
	}
	fclose(f);
}
