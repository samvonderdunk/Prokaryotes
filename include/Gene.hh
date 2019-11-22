#ifndef GeneHeader
#define GeneHeader

#include "Bead.hh"

class Gene : public Bead {
 public:
   int threshold; //Note that a threshold actually represents the lowest value for which a gene becomes active, i.e. if threshold=1, then the actual threshold lies at 0.5 and the gene will be expressed if the activity sums up to 1 or more.
   int activity;  //TF behaviour depends on the TF, i.e. the gene; not the TFBS.
   bool binding_domain[binding_length];
   int expression;

  Gene();//constructor
  Gene(int t, int th, int act, bool bd_dom[], int expr);
  explicit Gene(const Gene &gene);//copy constructor
  ~Gene();//destructor

  virtual Bead* Clone() const;
  void RandomGene();
};

#endif
