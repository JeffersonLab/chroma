/*! File: stag_scalars_s.cc 
 *
 * The routines in this file compute all 16 staggered scalars.
 * 
 * BEWARE: These routines ASSUME that the pion propagators have been
 * calculated in the Coulomb (or other spatially fixed) gauge.
 * It is not gauge invariant. It could be made to be so by adding some
 * parallel transport, however folklore claims that increases noise
 *
 * YOU HAVE BEEN WARNED.
 */

#include "meas/hadron/stag_propShift_s.h"
#include "meas/hadron/stag_scalars_s.h"
#include "util/gauge/stag_phases_s.h"


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
using namespace StagPhases;

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
  UnorderedSet timeslice;
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
  corr_fn_s = - alpha(1)*beta(0)*trace(adj(quark_props[0])*quark_props[0]);

  // Slice Sum
  corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);

  sca_index++;


  // Array to describe shifts in cube
  multi1d<int> delta(Nd);

  // One link spatial scalars.
  for(mu=0; mu<Nd-1; mu++) {
 
    delta = 0;
    delta[mu] = 1;
      
    corr_fn_s =  alpha(mu+1)*trace(shift_deltaProp(delta,quark_props[0])
                             *adj(quark_props[ deltaToPropIndex(delta) ]));

    corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
    sca_index++;
  }
    
  // One link temporal
  corr_fn_s = -  beta(0)*alpha(1)*alpha(3)*trace(adj(quark_props[0])*quark_props[0]);
  corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);

  sca_index++;

  // Two link spatial
  for(mu=0; mu<Nd-1; mu++) { 
    for(nu=mu+1; nu <Nd-1; nu++) { 
      delta = 0;
      delta[mu] = 1;
      delta[nu] = 1;

      corr_fn_s = beta(mu)* alpha(nu+1)
                              *trace(adj(shift_deltaProp(delta,quark_props[0]))
                                     *quark_props[ deltaToPropIndex(delta) ]);

      corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
      sca_index++;
    }
  }

  // Two link temporal
  for(mu=0; mu<Nd-1; mu++) { 
      delta = 0;
      delta[mu] = 1;
    
      corr_fn_s = - beta(mu)*  beta(2)
                              *trace(adj(shift_deltaProp(delta,quark_props[0]))
                                     *quark_props[ deltaToPropIndex(delta) ]);

    
      corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
      sca_index++;
  }

  // Three link spatial
  for(mu=0; mu<Nd-1; mu++) { 
    for(nu=mu+1; nu <Nd-1; nu++) { 
      for(rho=nu+1; rho < Nd-1; rho++) { 

	delta = 0;
	delta[mu] = 1;
	delta[nu] = 1;
	delta[rho] = 1;

        corr_fn_s = - alpha(mu+1) * alpha(nu+1)* alpha(rho+1)
          *trace(adj(shift_deltaProp(delta,quark_props[0]))
                 *quark_props[ deltaToPropIndex(delta) ]);

        corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
	sca_index++;
      }
    }
  }
  
  // Three link temporal 
  for(mu=0; mu<Nd-1; mu++) { 
    for(nu=mu+1; nu <Nd-1; nu++) { 

	delta = 0;
	delta[mu] = 1;
	delta[nu] = 1;

	corr_fn_s =  beta(mu)* beta(nu)*  beta(2)
          *trace(adj(shift_deltaProp(delta,quark_props[0]))
                 *quark_props[ deltaToPropIndex(delta) ]);

	corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
	sca_index++;
    }
  }

  // Four link temporal
  delta = 0;
  delta[0] = delta[1] = delta[2] = 1;

  corr_fn_s = beta(0)* beta(1)
    *trace(adj(shift_deltaProp(delta, quark_props[0]))
           *quark_props[ deltaToPropIndex(delta) ] );

  corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
  sca_index++;

  if( sca_index != NUM_STAG_PIONS) { 
    QDPIO::cerr << "Panic! Panic! Something has gone horribly wrong" << endl;
    QDP_abort(1);
  }
}
