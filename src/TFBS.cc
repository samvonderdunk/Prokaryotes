#include "TFBS.hh"

TFBS::TFBS() : Bead() {	// TFBS is derived class of Bead
  type=1;
  activity=0;
  for(int i=0; i<binding_length; i++) binding_site[i] = 0;
  ClaimVector = new vector<double>();
}

TFBS::TFBS(int t, int act, bool bd_dom[]) : Bead() {
  type=t;
  activity=act;
  for(int i=0; i<binding_length; i++) binding_site[i] = bd_dom[i];
  ClaimVector = new vector<double>();
}

TFBS::TFBS(const TFBS &tfbs) : Bead(tfbs) {
  type=tfbs.type;
  activity=tfbs.activity;
  for(int i=0; i<binding_length; i++) binding_site[i] = tfbs.binding_site[i];
  ClaimVector = new vector<double>();
  vector<double>::iterator iv;
  iv = tfbs.ClaimVector->begin();
  while (iv != tfbs.ClaimVector->end())
  {
    ClaimVector->push_back(*iv);
    iv++;
  }
}

TFBS::~TFBS() {
  delete ClaimVector;
  ClaimVector=NULL;
}

Bead* TFBS::Clone() const {
  return new TFBS(*this);
}

void TFBS::RandomTFBS()
{
  type = 1;
  activity = (uniform()>0.5) ? -1 : 1;
  for (int k=0; k<binding_length; k++)	binding_site[k] = (uniform()>0.5) ? true : false;
}
