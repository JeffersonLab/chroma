// $Id: scalar_loops_s.cc,v 1.6 2005-06-28 16:27:21 mcneile Exp $

#include "chromabase.h"
#include"scalar_loops_s.h"
#include "util/gauge/stag_phases_s.h"

namespace Chroma {

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
  corr += corr_fn[isample] ;

}



void non_local_scalar_loop::compute(LatticeStaggeredFermion & q_source, 
				LatticeStaggeredFermion & psi, int isample)
{
  UnorderedSet timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  LatticeStaggeredFermion psi_sca1 ;

  // Array to describe shifts in cube
  multi1d<int> delta(Nd);
  delta = 0;
  delta[3] = 1;
  psi_sca1 = shift_deltaProp(delta, psi); 
  LatticeComplex TrG_s1 ; 

  using namespace StagPhases;
  TrG_s1 = alpha(3)*localInnerProduct(q_source, psi_sca1);

  corr_fn[isample] = sumMulti(TrG_s1, timeslice);
  corr += corr_fn[isample] ;

}

}  // end namespace Chroma
