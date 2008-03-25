/*! File: stag_scalars_s.cc 
 *
 * The routines in this file compute eight  1xTASTE and eight gamma3xTASTE
 * (gamma_txTASTE)
 * staggered scalars.
 * The latter may be exotic and the former may be redundant. Run it if you 
 * want. Your call.
 *
 * comments added to elucidate what operators are being computed
 * for consistency we keep the same crazy indexing, and for now the same 
 * XML tags. This whole file might be redundant.
 *
 */

#include "chromabase.h"
#include "meas/hadron/stag_propShift_s.h"
#include "meas/hadron/stag_scalars_s.h"
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

  /*! staggeredScalars
   *
   * This routine computes all 16 staggered scalars.
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
   *       scalar_corr_fn     -- The 2d array of scalar correlation functions
   *                           (16, Nt)  (Write)
   *
   *       j_decay          -- The time direction (has to be Nd-1 for now)
   *                           (Read)
   */
 

// THis brings the staggered phases alpha and beta into the namespace

  void 
  staggered_scalars::compute(
    multi1d<LatticeStaggeredPropagator>& quark_props,
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
      QDPIO::cerr << "staggeredScalars: input quark props has the wrong number of elements. It should be 8 but is " << quark_props.size() << endl;
      QDP_abort(1);
    };

    // Also for now we want j_decay to be 3.
    switch( j_decay ) { 
    case 3:
      break;
    
    default:
      QDPIO::cerr << "staggeredScalars: j_decay must be 3 for just now. It is " << j_decay << endl;
      QDP_abort(1);
    };

    // Get the lattice size.
    const multi1d<int>& latt_size = Layout::lattSize();
  
    // resize output array appropriately
    corr_fn.resize(NUM_STAG_PIONS, latt_size[Nd-1]);

    // Correlation functions before spatial sum
    LatticeComplex corr_fn_s;

    // Machinery to do timeslice sums with 
    Set timeslice;
    timeslice.make(TimeSliceFunc(Nd-1));

    // Phases
    //multi1d<LatticeInteger> alpha(Nd); // KS Phases
    //multi1d<LatticeInteger> beta(Nd);  // Auxiliary phases for this work

    // Get the phases -- now done elsewhere
    // mesPhasFollana(alpha, beta);
  
    // Counters/Indices
    int sca_index = 0;
    int i;
    int mu, nu, rho;  

    // Taste singlet scalar (connected correlator)
    //  1x1
    corr_fn_s = - StagPhases::alpha(1)*StagPhases::beta(0)*trace(adj(quark_props[0])*quark_props[0]);

    // Slice Sum
    corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);

    sca_index++;


    // Array to describe shifts in cube
    multi1d<int> delta(Nd);

    // One link spatial scalars.
    // 1xgamma0, 1xgamma1, 1xgamma2

    for(mu=0; mu<Nd-1; mu++) {
 
      delta = 0;
      delta[mu] = 1;
      
      corr_fn_s =  StagPhases::alpha(mu+1)*trace(shift_deltaProp(delta,quark_props[0])
				     *adj(quark_props[ deltaToPropIndex(delta) ]));

      corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
      sca_index++;
    }
    
    // zero link gamma3 operator
    // gamma3xgamma3

    corr_fn_s = -  StagPhases::beta(0)*StagPhases::alpha(1)*StagPhases::alpha(3)*trace(adj(quark_props[0])*quark_props[0]);
    corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);

    sca_index++;

    // Two link spatial
    // 1xgamma0gamma1, 1xgamma0gamm2, 1xgamma1gamma2

    for(mu=0; mu<Nd-1; mu++) { 
      for(nu=mu+1; nu <Nd-1; nu++) { 
	delta = 0;
	delta[mu] = 1;
	delta[nu] = 1;

	corr_fn_s = StagPhases::beta(mu)* StagPhases::alpha(nu+1)
	  *trace(adj(shift_deltaProp(delta,quark_props[0]))
		 *quark_props[ deltaToPropIndex(delta) ]);

	corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
	sca_index++;
      }
    }

    // one link gamma3 operators
    // gamma3xgamma0gamma3, gamma3xgamma1gamma3, gamma3xgamma2gamma3
    for(mu=0; mu<Nd-1; mu++) { 
      delta = 0;
      delta[mu] = 1;
    
      corr_fn_s = - StagPhases::beta(mu)*  StagPhases::beta(2)
	*trace(adj(shift_deltaProp(delta,quark_props[0]))
	       *quark_props[ deltaToPropIndex(delta) ]);

    
      corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
      sca_index++;
    }

    // Three link spatial scalar
    // despite the loops, there can be only one
    // 1xgamma3gamma5

    for(mu=0; mu<Nd-1; mu++) { 
      for(nu=mu+1; nu <Nd-1; nu++) { 
	for(rho=nu+1; rho < Nd-1; rho++) { 

	  delta = 0;
	  delta[mu] = 1;
	  delta[nu] = 1;
	  delta[rho] = 1;

	  corr_fn_s = - StagPhases::alpha(mu+1) * StagPhases::alpha(nu+1)* StagPhases::alpha(rho+1)
	    *trace(adj(shift_deltaProp(delta,quark_props[0]))
		   *quark_props[ deltaToPropIndex(delta) ]);

	  corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
	  sca_index++;
	}
      }
    }
  
    // two-link gamma3 operators
    // gamma3xgamma2gamma5, gamma3xgamma1gamma5, gamma3xgamma0gamma5
    // 
    for(mu=0; mu<Nd-1; mu++) { 
      for(nu=mu+1; nu <Nd-1; nu++) { 

	delta = 0;
	delta[mu] = 1;
	delta[nu] = 1;

	corr_fn_s =  StagPhases::beta(mu)* StagPhases::beta(nu)*  StagPhases::beta(2)
          *trace(adj(shift_deltaProp(delta,quark_props[0]))
                 *quark_props[ deltaToPropIndex(delta) ]);

	corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
	sca_index++;
      }
    }

    // three-link gamma3 operator
    // gamma3xgamma5

    delta = 0;
    delta[0] = delta[1] = delta[2] = 1;

    corr_fn_s = StagPhases::beta(0)* StagPhases::beta(1)
      *trace(adj(shift_deltaProp(delta, quark_props[0]))
	     *quark_props[ deltaToPropIndex(delta) ] );

    corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
    sca_index++;

    if( sca_index != NUM_STAG_PIONS) { 
      QDPIO::cerr << "Panic! Panic! Something has gone horribly wrong" << endl;
      QDP_abort(1);
    }
  }

}  // end namespace Chroma
