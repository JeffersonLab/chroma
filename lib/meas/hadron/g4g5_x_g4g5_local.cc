/*! File: pion_local_s.cc
 *
 * The routines in this file compute the staggered pions
 *  (gamma_4 gamma_5 cross gamma_4 gamma_5)
 *
 *
 */

#include "meas/hadron/g4g5_x_g4g5_local.h"
#include "util/gauge/stag_phases_s.h"

namespace Chroma {

// I cant forward declare this for some reason
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

/*! compute method for the g4g5_x_g4g5_local_meson class
 *
 * This routine computes all 16 staggered pions.
 * 
 * Caveats: i) It assumes that the propagators you give 
 *             have been computed in some spatially fixed gauge
 *             eg the Coulomb Gauge. 
 *
 *
 *         iii) The assumption is that you are working in 4d 
 * 
 *  Parameters: 
 * 
 *       quark_props      -- The array of input propagators (Read)
 *       pion_corr_fn     -- The 2d array of pion correlation functions
 *                           (16, Nt)  (Write)
 *
 *       j_decay          -- The time direction (has to be Nd-1 for now)
 *                           (Read)
 */
 

// THis brings the staggered phases alpha and beta into the namespace
  using namespace StagPhases;

  void 
  g4g5_x_g4g5_local_meson::compute(
				LatticeStaggeredPropagator& quark_prop_A,
				LatticeStaggeredPropagator& quark_prop_B,
				int j_decay)
  {

    // Paranoid Checks

    if( Nd != 4 ) { 
      QDPIO::cerr << "The no of dimensions should be 4 for now. It is: " 
		  << Nd << endl;
      QDP_abort(1);
    }

    // Also for now we want j_decay to be 3.
    switch( j_decay ) { 
    case 3:
      break;
    
    default:
      QDPIO::cerr << "pions_s: j_decay must be 3 for just now. It is " << j_decay << endl;
      QDP_abort(1);
    };

    // Get the lattice size.
    const multi1d<int>& latt_size = Layout::lattSize();
  

    // Correlation functions before spatial sum
    LatticeComplex latt_corr_fn;

    // Machinery to do timeslice sums with 
    Set timeslice;
    timeslice.make(TimeSliceFunc(Nd-1));

    // Counters/Indices
    int pion_index = 0;
    int i;
    int mu, nu, rho;  

    // Goldstone Pion
    latt_corr_fn = -  alpha(Nd-1)*trace(adj(quark_prop_A)*quark_prop_B);

    // Slice Sum
    corr_fn[ pion_index ] = sumMulti(latt_corr_fn, timeslice);
    tag_names[pion_index] = "gamma4gamma5_CROSS_gamma4gamma5" ; 

    pion_index++;

    if( pion_index != no_pions ) { 
      QDPIO::cerr << "Panic! Panic! Something has gone horribly wrong" << endl;
      QDP_abort(1);
    }
  }

}  // end namespace Chroma

