#include "meas/hadron/scalar_local_s.h"
#include "util/gauge/stag_phases_s.h"

namespace Chroma {

  class TimeSliceFunc : public SetFunc {
  public:
    TimeSliceFunc(int dir) : dir_decay(dir) {}
    int operator() (const multi1d<int>& coordinate) const { return coordinate[dir_decay]; }
    int numSubsets() const { return Layout::lattSize()[dir_decay]; }

    int dir_decay;

  private:
    TimeSliceFunc() {}
  };

  using namespace StagPhases;

  void staggered_local_scalar::compute(LatticeStaggeredPropagator& quark_prop_A,
				       LatticeStaggeredPropagator& quark_prop_B, int j_decay) 
  {
    if (Nd != 4 ) {
      QDPIO::cerr << "The number of dimensions should be 4 (for now) - it is " << Nd << "\n";
      QDP_abort(1);
    }

    switch(j_decay) {
    case 3:
      break;

    default:
      QDPIO::cerr << "scalar_s: j_decay must be 3 for now, it is " << j_decay << "\n";
      QDP_abort(1);
    }

    const multi1d<int>& latt_size = Layout::lattSize();

    LatticeComplex latt_corr_fn;

    Set timeslice;
    timeslice.make(TimeSliceFunc(Nd-1));

    int idx = 0;
    int i, mu, nu, rho;


    latt_corr_fn = - alpha(1)*beta(0)*trace(adj(quark_prop_A)*quark_prop_B);

    corr_fn[idx] = sumMulti(latt_corr_fn, timeslice);
    tag_names[idx] = "one_CROSS_one";

    idx++;

    if ( idx != no_pions ) {
      QDPIO::cerr << "Panic! Something has gone horribly wrong!\n";
      QDP_abort(1);
    }
  }


    
} //end namespace
