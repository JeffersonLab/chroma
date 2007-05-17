// $Id: pseudoscalar_loops_s.cc,v 3.3 2007-05-17 15:29:23 egregory Exp $
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
  // Array to describe shifts in cube
  multi1d<int> delta(Nd);

  Set timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  LatticeStaggeredFermion psi_eta4 ;

  delta = 0;
  delta[0] = delta[1] = delta[2] = delta[3] = 1;
 
  psi_eta4 = shift_deltaProp(delta, psi);

  LatticeComplex TrG_eta4 ; 
  using namespace StagPhases;

  TrG_eta4 = -alpha(0)*alpha(1)*alpha(2)*alpha(3)*localInnerProduct(q_source, psi_eta4);

  corr_fn[isample] = sumMulti(TrG_eta4, timeslice);
  corr += corr_fn[isample] ;
}



void threelink_pseudoscalar_loop::compute(
					  LatticeStaggeredFermion & q_source, 
					  LatticeStaggeredFermion & psi, 
					  int isample)
{
  // Array to describe shifts in cube
  multi1d<int> delta(Nd);

  Set timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  LatticeStaggeredFermion  psi_eta3 ;

  delta = 0;
  delta[0] = delta[1] = delta[2]  = 1;
 
  psi_eta3 = shift_deltaProp(delta, psi);

  LatticeComplex TrG_eta3 ; 

  using namespace StagPhases;
  TrG_eta3 = alpha(0)*alpha(1)*alpha(2)*localInnerProduct(q_source, psi_eta3);

  corr_fn[isample] = sumMulti(TrG_eta3, timeslice);
  corr += corr_fn[isample] ;

}

void fourlink_pseudoscalar_kilcup_loop::compute(LatticeStaggeredFermion & psi,
						int isample, Real mass){

  // Array to describe shifts in cube
  multi1d<int> delta(Nd);

  Set timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  LatticeStaggeredFermion psi_eta4 ;

  delta = 0;
  delta[0] = delta[1] = delta[2] = delta[3] = 1;
 
  psi_eta4 = shift_deltaProp(delta, psi);

  LatticeComplex TrG_eta4 ; 
  using namespace StagPhases;

  TrG_eta4 = -2*mass*alpha(0)*alpha(1)*alpha(2)*alpha(3)
                    *localInnerProduct(psi_eta4, psi);

  corr_fn[isample] = sumMulti(TrG_eta4, timeslice);
  corr += corr_fn[isample] ;
}

void fourlink_pseudoscalar_kilcup_loop::compute(
						LatticeStaggeredFermion & q_source,
				LatticeStaggeredFermion & psi, int isample){}


void zerolink_pseudoscalar_loop::compute(
					  LatticeStaggeredFermion & q_source, 
					  LatticeStaggeredFermion & psi, 
					  int isample)
{
  // Array to describe shifts in cube
  multi1d<int> delta(Nd);

  Set timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  LatticeStaggeredFermion  psi_eta0 ;

  delta = 0;
  //  delta[0] = delta[1] = delta[2]  = 1;
 
  //  psi_eta3 = shift_deltaProp(delta, psi);
  psi_eta0 = psi;

  LatticeComplex TrG_eta0 ; 

  using namespace StagPhases;
  TrG_eta0 = alpha(0)*beta(1)*localInnerProduct(q_source, psi_eta0);

  corr_fn[isample] = sumMulti(TrG_eta0, timeslice);
  corr += corr_fn[isample] ;

}

void fourlink_pseudoscalar_loop_fuzz::compute(LatticeStaggeredFermion & q_source, 
				LatticeStaggeredFermion & psi_fuzz, int isample)
{
  // Array to describe shifts in cube
  multi1d<int> delta(Nd);

  Set timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  LatticeStaggeredFermion psi_eta4 ;

  delta = 0;
  delta[0] = delta[1] = delta[2] = delta[3] = 1;
 
  psi_eta4 = shift_deltaProp(delta, psi_fuzz);

  LatticeComplex TrG_eta4 ; 
  using namespace StagPhases;

  TrG_eta4 = -alpha(0)*alpha(1)*alpha(2)*alpha(3)*localInnerProduct(q_source, psi_eta4);

  corr_fn[isample] = sumMulti(TrG_eta4, timeslice);
  corr += corr_fn[isample] ;
}



void threelink_pseudoscalar_loop_fuzz::compute(
					  LatticeStaggeredFermion & q_source, 
					  LatticeStaggeredFermion & psi_fuzz, 
					  int isample)
{
  // Array to describe shifts in cube
  multi1d<int> delta(Nd);

  Set timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  LatticeStaggeredFermion  psi_eta3 ;

  delta = 0;
  delta[0] = delta[1] = delta[2]  = 1;
 
  psi_eta3 = shift_deltaProp(delta, psi_fuzz);

  LatticeComplex TrG_eta3 ; 

  using namespace StagPhases;
  TrG_eta3 = alpha(0)*alpha(1)*alpha(2)*localInnerProduct(q_source, psi_eta3);

  corr_fn[isample] = sumMulti(TrG_eta3, timeslice);
  corr += corr_fn[isample] ;

}

void fourlink_pseudoscalar_kilcup_loop_fuzz::compute(LatticeStaggeredFermion & psi_fuzz,
					       LatticeStaggeredFermion & psi,
						int isample, Real mass){

  // Array to describe shifts in cube
  multi1d<int> delta(Nd);

  Set timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  LatticeStaggeredFermion psi_eta4 ;

  delta = 0;
  delta[0] = delta[1] = delta[2] = delta[3] = 1;
 
  psi_eta4 = shift_deltaProp(delta, psi);

  LatticeComplex TrG_eta4 ; 
  using namespace StagPhases;

  TrG_eta4 = -2*mass*alpha(0)*alpha(1)*alpha(2)*alpha(3)
                    *localInnerProduct(psi_eta4, psi_fuzz);

  corr_fn[isample] = sumMulti(TrG_eta4, timeslice);
  corr += corr_fn[isample] ;
}

void fourlink_pseudoscalar_kilcup_loop_fuzz::compute(
						LatticeStaggeredFermion & q_source,
				LatticeStaggeredFermion & psi_fuzz, int isample){}


}  // end namespace Chroma
