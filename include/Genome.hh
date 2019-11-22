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

  // vector<bool>* MutationList;  //Checks if every bead has had the chance to be mutated.
  // typedef vector<bool>::iterator mut_iter;
  // int deletion_length;

  vector<int>* GeneStates;
  vector<int>* GeneTypes;
  typedef vector<int>::iterator gene_iter;

  int g_length;
  int gnr_genes;
  int pos_fork, pos_anti_ori;  //position of the fork moving backwards through the parent genome and position of the origin of replication equal to the starting position of the fork.
  bool mutant_genome;

  Genome();
  ~Genome();

  void CopyPartOfGenome(iter begin, iter end);
  void CopyPartOfGenomeToTemplate(iter begin, iter end, list<Bead*>* template_beadlist);
  void CloneGenome(const Genome* G_template);

  void DevelopChildrenGenomes(Genome* G_replicated);
  void SplitGenome(Genome* G_replicated);

  void ReplicateGenomeStep(int env);
  iter GeneMutate(iter ii, int* pdel_len);
  iter TFBSMutate(iter ii, int* pdel_len);

  iter GeneDuplication(iter ii, int* pdup_len);
  iter GeneInnovation();
  iter GeneDeletion(iter ii, int* pdel_len);
  iter GeneShuffle(iter ii);

  iter TFBSDuplication(iter ii);
  void TFBSInnovation();
  iter TFBSDeletion(iter ii);
  iter TFBSShuffle(iter ii);

  iter FindFirstTFBSInFrontOfGene(iter ii) const;
  iter FindRandomGenePosition() const;

  void PotentialTypeChange(iter ii);
  bool CheckSameGeneTypes(iter ii, iter jj);
  bool CheckSameGeneTypes(Gene* gene_ii, Gene* gene_jj);  //Some redundancy (forgot the coding jargon), allowing me to also check the types of two gene pointers.
  int CountTypeAbundance(int type);

  void ReadInitialGenome();
  void ReadBeadsFromString(string genome);
  void ReadInitialGeneStates();
  void InitialiseRandomGenome();

  int MatchNextState(int state);
  void UpdateGeneStates();
  void SetClaimVectors();

  iter MatchGeneToTFBS(iter i_tfbs);
  int MatchBitStrings(Bead* b_tfbs, Bead* b_gene);

  void GenomeToNetwork(double** Network);

  bool IsGene(Bead* bead) const;
  bool IsTFBS(Bead* bead) const;
  void DecrementExpressionOfType(iter ii);
  void IncrementExpressionOfType(iter ii);
  int FindIndexOfType(int type);

  string PrintContent(list<Bead*>* chromosome, bool terminal, bool only_parent);
  string PrintGeneStateContent(bool sorted);
  string PrintGeneTypeContent();
};
#endif
