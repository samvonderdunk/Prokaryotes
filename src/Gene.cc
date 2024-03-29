#include "Gene.hh"
#include "Header.hh"

Gene::Gene() : Bead()
{
	int i;
  type = 0;
  threshold = 0;
  activity = 0;
	for(i=0; i<typeseq_length; i++)	typeseq[i] = 0;
  for(i=0; i<binding_length; i++) binding_domain[i] = 0;
  expression = 0;
}

Gene::Gene(int t, int th, int act, bool tseq[], bool bd_dom[], int expr) : Bead()
{
	int i;
  type=t;
  threshold=th;
  activity=act;
	for(i=0; i<typeseq_length; i++) typeseq[i] = tseq[i];
  for(i=0; i<binding_length; i++) binding_domain[i] = bd_dom[i];
  expression=expr;
}

Gene::Gene(const Gene &gene) : Bead(gene)    //See below.. indeed this function copies from the object gene (class Gene; not a pointer or anything!), which you can tell from the fact that members of gene are accessed by dots, i.e. gene.type
{
	int i;
  type=gene.type;
  threshold=gene.threshold;
  activity=gene.activity;
	for(i=0; i<typeseq_length; i++) typeseq[i] = gene.typeseq[i];
  for(i=0; i<binding_length; i++) binding_domain[i] = gene.binding_domain[i];
  expression=gene.expression;
}

Gene::~Gene() {
}

Bead* Gene::Clone() const    //This function returns a pointer to a bead (a gene in this case). It creates a pointer to a Gene via p_Gene=new Gene(); this is how you can initialise a new pointer to the type Gene. However, it does not call the constructor, because the brackets are not empty. Instead it calls the copy constructor (see above), which by rule requires the reference to a Gene class object (the ampersand does not indicate the "address of gene" but "reference to gene"?). It copies from itself: this is a pointer to the current object, so by dereferencing to *this, we are passing the object itself to the copy constructor. See above.
{
  return new Gene(*this);   //Note that you have to return a pointer to a Bead type, because Bead is only an abstract type. With the pointer to Bead type that this function makes/returns, you in this case also get the Gene type members the current object is a Gene and so *this is also a Gene. However doing it in this manner means that the function CopyPartOfGenome in Genome.cc can apply to any bead (any class derived from the Bead type). So you generally use pointers to the base class (which drag along their derived class).
}

void Gene::RandomGene()
{
	int k;
	type = 0;
  threshold = (int)(uniform()*(2*WeightRange+1) - WeightRange);	//Threshold between -3 and 3 (including these borders).
  activity = (int)(uniform()*(2*WeightRange+1) - WeightRange);
	for (k=0; k<typeseq_length; k++)	typeseq[k] = (uniform()>0.5) ? true : false;
  for (k=0; k<binding_length; k++)	binding_domain[k] = (uniform()>0.5) ? true : false;
  expression = 0;
}
