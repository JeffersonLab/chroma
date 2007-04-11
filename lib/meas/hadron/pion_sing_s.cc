/*! File:  pion_sing_s.cc   
 *
 * The routines in this file computes the four-link flavor-singlet pion
 * (This routine has not yet been CHECKED !!!!!)
 * 
 * BEWARE: These routines ASSUME that the quark propagators have been
 * calculated in the Coulomb (or other spatially fixed) gauge.
 * It is not gauge invariant. It could be made to be so by adding some
 * parallel transport, however folklore claims that increases noise
 *
 * YOU HAVE BEEN WARNED.
 */

#include "chromabase.h"
#include "pion_sing_s.h"
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

  /*! pion_singlet
   *
   * This routine computes 1 staggered four-link flavor-singlet pion.
   * 
   * Caveats: i) It assumes that the propagators you give 
   *             have been computed in some spatially fixed gauge
   *             eg the Coulomb Gauge. 
   *
   *          ii) This means that there is only 
   *             1 propagators corresponding to 
   *                  prop_index = 0,   hypercube_coord (1,1,1,1)
   *
   *         iii) The assumption is that you are working in 4d 
   * 
   *  Parameters: 
   * 
   *       quark_props      -- The array of input propagators (Read)
   *
   *       j_decay          -- The time direction (has to be Nd-1 for now)
   *                           (Read)
   */
 

  // THis brings the staggered phases alpha and beta into the namespace 
  using namespace StagPhases;

  void 
  staggered_pion_singlet::compute(
    LatticeStaggeredPropagator local_quark_prop,
    LatticeStaggeredPropagator four_shift_quark_prop,
    int j_decay)
  {

    //non-fuzzed version

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
      QDPIO::cerr << "pion_sing_s: j_decay must be 3 for just now. It is " << j_decay << endl;
      QDP_abort(1);
    };

    // Get the lattice size.
    const multi1d<int>& latt_size = Layout::lattSize();
  
    // resize output array appropriately
    corr_fn.resize(no_pion_sings, latt_size[Nd-1]);

    // Correlation functions before spatial sum
    LatticeComplex corr_fn_s;

    // Machinery to do timeslice sums with 
    Set timeslice;
    timeslice.make(TimeSliceFunc(Nd-1));

    // Phases
    // Counters/Indices
    int pion_index = 0;
    int i;
    int mu, nu, rho;  

    // 
    // Array to describe shifts in cube
    multi1d<int> delta(Nd);


    //---------------------------------------
    // four-link taste singlet pion

    tag_names[pion_index]     = "gamma5_CROSS_one" ; 

    delta = 0;
    delta[0] = delta[1] = delta[2] = delta[3] = 1;
    corr_fn_s = beta(0)* beta(1) * beta(2) * beta(3)
      *trace(adj(shift_deltaProp(delta, local_quark_prop ))
	     *four_shift_quark_prop);
  
    corr_fn[ pion_index ] = sumMulti(corr_fn_s, timeslice);
    pion_index++;

    if( pion_index !=  no_pion_sings ) { 
      QDPIO::cerr << "Panic! Panic! Something has gone horribly wrong" << endl;
      QDP_abort(1);
    }



  }// end non-fuzzed sink version

  /**********************************************************************/

}  // end namespace Chroma
