#ifndef TFBSHeader
#define TFBSHeader

#include "Bead.hh"

class TFBS : public Bead {
 public:
  bool binding_site[binding_length]; 	// Ranges from -weightrange to weightrange
  int activity; //This could be an element of the basal bead type, but than recoding might be more difficult in the future.
  vector<double>* ClaimVector;   //Each element corresponds to an element in GeneStates and in GeneTypes.

  typedef vector<double>::iterator claim_iter;

  TFBS();//constructor
  TFBS(int t, int act, bool bd_dom[]);
  explicit TFBS(const TFBS &tfbs);//copy constructor
  virtual ~TFBS();//destructor

  virtual Bead* Clone() const;
  void RandomTFBS();
};

#endif
