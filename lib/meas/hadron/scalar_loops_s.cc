
#include "chroma.h"
#include"scalar_loops_s.h"
#include "util/gauge/stag_phases_s.h"

// Standard Time Slicery
class TimeSliceFunc : public SetFunc
{
public:
  TimeSliceFunc(int dir): dir_decay(dir) {}
 
  int operator() (const multi1d<int>& coordinate) const {return coordinate[dir_decay];}
  int numSubsets() const {return Layout::lattSize()[dir_decay];}
 
  int dir_decay;
 
private:
  TimeSliceFunc() {}  // hide default constructor
};


void local_scalar_loop::compute(LatticeStaggeredFermion & q_source, 
				LatticeStaggeredFermion & psi, int isample)
{
  UnorderedSet timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  LatticeComplex TrG_s0 ; 
  TrG_s0 = localInnerProduct(q_source, psi);

  corr_fn[isample] = sumMulti(TrG_s0, timeslice);

}
