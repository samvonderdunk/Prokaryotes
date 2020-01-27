#include "House.hh"
#include "Header.hh"

House::House() : Bead()
{
  type = 1;
}

House::~House() {
}

House::House(const House &house) : Bead(house)
{
  type=house.type;
}

Bead* House::Clone() const
{
  return new House(*this);
}
