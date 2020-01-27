#ifndef HouseHeader
#define HouseHeader

#include "Bead.hh"

//Class House defines household genes that are impossible or deleterious to add/remove.

class House : public Bead {
 public:

  House();//constructor
  ~House();//destructor
  explicit House(const House &house);
  virtual Bead* Clone() const;

};

#endif
