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

  vector<bool>* MutationList;  //Checks if every bead has had the chance to be mutated.
  typedef vector<bool>::iterator mut_iter;
  int deletion_length;

  vector<int>* GeneStates;
  vector<int>* GeneTypes;
  typedef vector<int>::iterator gene_iter;

  int g_length;
  int gnr_genes;
  int pos_fork, pos_anti_ori;  //position of the fork moving backwards through the parent genome and position of the origin of replication equal to the starting position of the fork.

  Genome();
  ~Genome();

  void CopyPartOfGenome(iter begin, iter end);
  void CopyPartOfGenomeToTemplate(iter begin, iter end, list<Bead*>* template_beadlist);
  void CloneGenome(const Genome* G_template);

  void RemoveGenomeInParent(iter begin, iter end);
  void MoveGenomeToChild(iter begin, iter end);
  void SplitGenome(Genome* G_replicated);

  bool ReplicateGenomeStep();
  bool GeneMutate(iter ii);
  bool TFBSMutate(iter ii);

  iter GeneDuplication(iter ii, int* pdup_len);
  void GeneDeletion(iter ii);
  iter TFBSDuplication(iter ii);
  void TFBSDeletion(iter ii);
  iter FindFirstTFBSInFrontOfGene(iter ii) const;
  iter FindRandomGenePosition() const;

  void PotentialTypeChange(iter ii);
  bool CheckSameGeneTypes(iter ii, iter jj);
  void LoseGeneType(int type);
  int CountTypeAbundance(int type);

  void ReadInitialGenome();
  void ReadBeadsFromString(string genome);
  void ReadInitialGeneStates();
  void InitialiseRandomGenome();

  void UpdateGeneStates();
  void SetClaimVectors();

  iter MatchGeneToTFBS(iter i_tfbs);
  double MatchBitStrings(Bead* b_tfbs, Bead* b_gene);

  void GenomeToNetwork(double** Network);

  bool IsGene(Bead* bead) const;
  bool IsTFBS(Bead* bead) const;
  void DecrementExpressionOfType(iter ii);
  void IncrementExpressionOfType(iter ii);
  int FindIndexOfType(int type);

  string PrintContent(list<Bead*>* chromosome, bool terminal, bool only_parent);
  string PrintGeneStateContent();
  string PrintGeneTypeContent();
};
#endif
