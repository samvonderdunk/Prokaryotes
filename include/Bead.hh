#ifndef BeadHeader
#define BeadHeader

#include "Header.hh"

class Bead {
 public:
  int type;

  Bead();
  explicit Bead(const Bead &cbb);
  virtual ~Bead();

  virtual Bead* Clone() const=0;
};

#endif
