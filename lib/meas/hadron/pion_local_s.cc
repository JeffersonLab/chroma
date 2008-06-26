/*! File: pion_local_s.cc
 *
 * The routines in this file compute the local staggered pions
 *  (gamma_5 cross gamma_5)
 *
 *
 */

#include "meas/hadron/pion_local_s.h"
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

/*! compute method for the staggered_local_pion class
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
  staggered_local_pion::compute(
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
    latt_corr_fn = trace(adj(quark_prop_A)*quark_prop_B);

    // Slice Sum
    corr_fn[ pion_index ] = sumMulti(latt_corr_fn, timeslice);
    tag_names[pion_index] = "gamma5_CROSS_gamma5" ; 

    pion_index++;

    if( pion_index != no_pions ) { 
      QDPIO::cerr << "Panic! Panic! Something has gone horribly wrong" << endl;
      QDP_abort(1);
    }
  }


  /**
       Include the momentum phase factors


  **/

  void 
  staggered_local_pion::compute_and_dump(
				LatticeStaggeredPropagator& quark_prop_A,
				LatticeStaggeredPropagator& quark_prop_B,
				int j_decay, int t0, 
				const SftMom& phases, XMLWriter& xml)
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
    latt_corr_fn = trace(adj(quark_prop_A)*quark_prop_B);

    multi2d<DComplex> hsum;
    hsum = phases.sft(latt_corr_fn);

    XMLArrayWriter xml_sink_mom(xml,phases.numMom());

    push(xml_sink_mom, "gamma5_CROSS_gamma5_mom");

    for (int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
    {
      push(xml_sink_mom);
      write(xml_sink_mom, "sink_mom_num", sink_mom_num);
      write(xml_sink_mom, "sink_mom", phases.numToMom(sink_mom_num));

      int tt_length = latt_size[Nd-1]  ; // DEBUG 
      multi1d<DComplex> mesprop(tt_length);
      for (int t=0; t < tt_length; ++t) 
      {
        int t_eff = (t - t0 + tt_length) % tt_length;
	mesprop[t_eff] = hsum[sink_mom_num][t];
      }

      write(xml_sink_mom, "mesprop", mesprop);
      pop(xml_sink_mom);

    } // end for(sink_mom_num)
 
    //    push(xml_sink_mom);



  }



}  // end namespace Chroma

