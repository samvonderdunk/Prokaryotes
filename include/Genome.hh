#ifndef GenomeHeader
#define GenomeHeader

#include "Header.hh"
#include "Bead.hh"
#include "Gene.hh"
#include "TFBS.hh"
#include <typeinfo>

class Genome {
 public:
  std::list<Bead*>* BeadList;
  typedef std::list<Bead*>::iterator iter;
  typedef std::list<Bead*>::reverse_iterator reviter;

  vector<int>* GeneStates;
  vector<int>* GeneTypes;
  typedef vector<int>::iterator gene_iter;

  int g_length;
  int gnr_genes;
  bool mutated_bitstring;
  bool mutated_type;
  bool mutated_new_tfbs;

  Genome();
  ~Genome();

  void CopyPartOfGenome(iter begin, iter end);
  void CopyPartOfGenomeToTemplate(iter begin, iter end, list<Bead*>* template_beadlist);
  void CloneGenome(const Genome* G_template);
  bool ReplicateGenome();

  void MutateGenome();
  iter GeneMutate(iter ii);
  void PotentialTypeChange(iter ii);
  iter TFBSMutate(iter ii);

  iter GeneDuplication(iter ii);
  iter GeneDeletion(iter ii);
  iter TFBSDuplication(iter ii);
  iter TFBSDeletion(iter ii);
  iter FindFirstTFBSInFrontOfGene(iter ii) const;
  iter FindRandomGenePosition() const;

  bool CheckSameGeneTypes(iter ii, iter jj);
  void LoseGeneType(int type);
  int CountTypeAbundance(int type);

  void ReadInitialGenome();
  void ReadInitialGeneStates();
  void InitialiseRandomGenome();

  void UpdateGeneStates();
  void SetClaimVectors();

  iter MatchGeneToTFBS(iter i_tfbs);
  double MatchBitStrings(Bead* b_tfbs, Bead* b_gene);
  int UpdateSingleGene(iter ii);

  bool IsGene(Bead* bead) const;
  bool IsTFBS(Bead* bead) const;
  string PrintContent(list<Bead*>* chromosome, bool terminal);
  string PrintGeneStateContent();
  string PrintGeneTypeContent();
};
#endif
