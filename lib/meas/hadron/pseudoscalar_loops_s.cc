// $Id: pseudoscalar_loops_s.cc,v 1.3 2005-01-14 18:42:36 edwards Exp $
#include "chromabase.h"
#include "pseudoscalar_loops_s.h"
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


void fourlink_pseudoscalar_loop::compute(LatticeStaggeredFermion & q_source, 
				LatticeStaggeredFermion & psi, int isample)
{
  UnorderedSet timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  LatticeStaggeredFermion psi_eta4 ;
  psi_eta4 = shift(shift(shift(shift(psi, FORWARD, 0), FORWARD, 1),
			 FORWARD, 2), FORWARD, 3);
  LatticeComplex TrG_eta4 ; 
  using namespace StagPhases;

  TrG_eta4 = -alpha(1)*alpha(2)*alpha(3)*localInnerProduct(q_source, psi_eta4);

  corr_fn[isample] = sumMulti(TrG_eta4, timeslice);
  corr += corr_fn[isample] ;
}



void threelink_pseudoscalar_loop::compute(
					  LatticeStaggeredFermion & q_source, 
					  LatticeStaggeredFermion & psi, 
					  int isample)
{
  UnorderedSet timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  LatticeStaggeredFermion  psi_eta3 ;
  psi_eta3 = shift(shift(shift(psi, FORWARD, 0), FORWARD, 1), FORWARD, 2);

  LatticeComplex TrG_eta3 ; 

  using namespace StagPhases;
  TrG_eta3 = alpha(3)*localInnerProduct(q_source, psi_eta3);

  corr_fn[isample] = sumMulti(TrG_eta3, timeslice);
  corr += corr_fn[isample] ;

}

}  // end namespace Chroma
