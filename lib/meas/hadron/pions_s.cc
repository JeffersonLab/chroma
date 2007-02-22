/*! File: pions_follana_s.cc 
 *
 * The routines in this file compute all 16 staggered pions.
 *
 * BEWARE: These routines ASSUME that the pion propagators have been
 * calculated in the Coulomb (or other spatially fixed) gauge.
 * It is not gauge invariant. It could be made to be so by adding some
 * parallel transport, however folklore claims that increases noise
 *
 * YOU HAVE BEEN WARNED.
 */

#include "meas/hadron/pions_s.h"
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

/*! staggeredPionsFollana
 *
 * This routine computes all 16 staggered pions.
 * 
 * Caveats: i) It assumes that the propagators you give 
 *             have been computed in some spatially fixed gauge
 *             eg the Coulomb Gauge. 
 *
 *          ii) This means that there is only 
 *             8 propagators corresponding to the 8 corners of the 
 *             spatial cube. The props come in an array whose single
 *             index maps lexicographically to the corners of the cube.
 *             ie:  prop_index = 0,   hypercube_coord (0,0,0,0)
 *                  prop_index = 1,   hypercube_coord (1,0,0,0)
 *                  prop_index = 2,   hypercube_coord (0,1,0,0)
 *
 * essentially prop_index = x + 2*y + 4*z
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
  staggered_pions::compute(multi1d<LatticeStaggeredPropagator>& quark_props,
			   //		      multi2d<DComplex>& pion_corr_fn,
			   int j_decay)
  {

    // Paranoid Checks

    if( Nd != 4 ) { 
      QDPIO::cerr << "The no of dimensions should be 4 for now. It is: " 
		  << Nd << endl;
      QDP_abort(1);
    }

    // Check for 8 props
    if( quark_props.size() != NUM_STAG_PROPS ) { 
      QDPIO::cerr << "pions_s: input quark props has the wrong number of elements. It should be 8 but is " << quark_props.size() << endl;
      QDP_abort(1);
    };

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
  
    // resize output array appropriately
    corr_fn.resize(NUM_STAG_PIONS, latt_size[Nd-1]);

    // Correlation functions before spatial sum
    LatticeComplex latt_corr_fn;

    // Machinery to do timeslice sums with 
    Set timeslice;
    timeslice.make(TimeSliceFunc(Nd-1));

    // Phases
    //multi1d<LatticeInteger> alpha(Nd); // KS Phases
    //multi1d<LatticeInteger> beta(Nd);  // Auxiliary phases for this work

    // Get the phases -- now done elsewhere
    // mesPhasFollana(alpha, beta);
  
    // Counters/Indices
    int pion_index = 0;
    int i;
    int mu, nu, rho;  

    // Goldstone Pion
    latt_corr_fn = trace(adj(quark_props[0])*quark_props[0]);

    // Slice Sum
    corr_fn[ pion_index ] = sumMulti(latt_corr_fn, timeslice);
    tag_names[pion_index] = "gamma5_CROSS_gamma5" ; 

    pion_index++;


    // Array to describe shifts in cube
    multi1d<int> delta(Nd);

    tag_names[pion_index]   = "gamma5_CROSS_gamma0_gamma5" ; 
    tag_names[pion_index+1] = "gamma5_CROSS_gamma1_gamma5" ; 
    tag_names[pion_index+2] = "gamma5_CROSS_gamma2_gamma5" ; 
    for(mu=0; mu<Nd-1; mu++) {
 
      delta = 0;
      delta[mu] = 1;
      
      latt_corr_fn =  beta(mu)*trace(shift_deltaProp(delta,quark_props[0])
				     *adj(quark_props[ deltaToPropIndex(delta) ]));
    
      corr_fn[ pion_index ] = sumMulti(latt_corr_fn, timeslice);
      pion_index++;
    }
    
    // ------------------------------
    tag_names[pion_index]   = "gamma3_gamma5_CROSS_gamma3_gamma5" ; 
    latt_corr_fn = -  alpha(Nd-1)*trace(adj(quark_props[0])*quark_props[0]);
    corr_fn[ pion_index ] = sumMulti(latt_corr_fn, timeslice);
    pion_index++;

    // -----------------------------
    tag_names[pion_index]     = "gamma5_CROSS_gamma2_gamma3" ; 
    tag_names[pion_index+1]   = "gamma5_CROSS_gamma1_gamma3" ; 
    tag_names[pion_index+2]   = "gamma5_CROSS_gamma0_gamma3" ; 

    for(mu=0; mu<Nd-1; mu++) { 
      for(nu=mu+1; nu <Nd-1; nu++) { 
	delta = 0;
	delta[mu] = 1;
	delta[nu] = 1;

	latt_corr_fn = - beta(mu)* beta(nu)
	  *trace(adj(shift_deltaProp(delta,quark_props[0]))
		 *quark_props[ deltaToPropIndex(delta) ]);
    
	corr_fn[ pion_index ] = sumMulti(latt_corr_fn, timeslice);
	pion_index++;
      }
    }

    // --------------------

    tag_names[pion_index]     = "gamma3_gamma5_CROSS_gamma1_gamma2" ; 
    tag_names[pion_index+1]   = "gamma3_gamma5_CROSS_gamma0_gamma2" ; 
    tag_names[pion_index+2]   = "gamma3_gamma5_CROSS_gamma0_gamma1" ; 

    for(mu=0; mu<Nd-1; mu++) { 
      delta = 0;
      delta[mu] = 1;
    
      latt_corr_fn = - beta(mu)*  alpha(Nd-1)
	*trace(adj(shift_deltaProp(delta,quark_props[0]))
	       *quark_props[ deltaToPropIndex(delta) ]);
    
      corr_fn[ pion_index ] = sumMulti(latt_corr_fn, timeslice);
      pion_index++;
    }

    // ---------------------------------
    tag_names[pion_index]     = "gamma5_CROSS_gamma3" ; 

    for(mu=0; mu<Nd-1; mu++) { 
      for(nu=mu+1; nu <Nd-1; nu++) { 
	for(rho=nu+1; rho < Nd-1; rho++) { 

	  delta = 0;
	  delta[mu] = 1;
	  delta[nu] = 1;
	  delta[rho] = 1;

	
	  latt_corr_fn = - beta(mu) * beta(nu)* beta(rho)
	    *trace(adj(shift_deltaProp(delta,quark_props[0]))
		   *quark_props[ deltaToPropIndex(delta) ]);
	
	  corr_fn[ pion_index ] = sumMulti(latt_corr_fn, timeslice);
	  pion_index++;
	}
      }
    }
  
    // ---------------------------------
    tag_names[pion_index]     = "gamma3_gamma5_CROSS_gamma2" ; 
    tag_names[pion_index+1]   = "gamma3_gamma5_CROSS_gamma1" ; 
    tag_names[pion_index+2]   = "gamma3_gamma5_CROSS_gamma0" ; 

    for(mu=0; mu<Nd-1; mu++) { 
      for(nu=mu+1; nu <Nd-1; nu++) { 

	delta = 0;
	delta[mu] = 1;
	delta[nu] = 1;

	latt_corr_fn =  beta(mu)* beta(nu)*  alpha(Nd-1)
	  *trace(adj(shift_deltaProp(delta,quark_props[0]))
		 *quark_props[ deltaToPropIndex(delta) ]);
	
	corr_fn[ pion_index ] = sumMulti(latt_corr_fn, timeslice);
	pion_index++;
      }
    }

    // ---------------------------------
    tag_names[pion_index]     = "gamma3_gamma5_CROSS_one" ; 

    delta = 0;
    delta[0] = delta[1] = delta[2] = 1;
    latt_corr_fn = - alpha(3)* beta(0)* beta(1) * beta(2)
      *trace(adj(shift_deltaProp(delta, quark_props[0]))
	     *quark_props[ deltaToPropIndex(delta) ] );
  
    corr_fn[ pion_index ] = sumMulti(latt_corr_fn, timeslice);
    pion_index++;

    if( pion_index != NUM_STAG_PIONS) { 
      QDPIO::cerr << "Panic! Panic! Something has gone horribly wrong" << endl;
      QDP_abort(1);
    }
  }

}  // end namespace Chroma

