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
#include "chroma.h"
// #include "mesphas_follana_s.h"
#include "pions_follana_s.h"
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

// Forward Declarations:
int 
deltaToPropIndex(multi1d<int>& delta);

LatticePropagator 
shiftDeltaProp(multi1d<int>& delta, const LatticePropagator& src);

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
staggeredPionsFollana(multi1d<LatticePropagator>& quark_props,
		      multi2d<DComplex>& pion_corr_fn,
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
    QDPIO::cerr << "staggeredPionsFollana: input quark props has the wrong number of elements. It should be 8 but is " << quark_props.size() << endl;
    QDP_abort(1);
  };

  // Also for now we want j_decay to be 3.
  switch( j_decay ) { 
  case 3:
    break;
    
  default:
    QDPIO::cerr << "staggeredPionsFollana: j_decay must be 3 for just now. It is " << j_decay << endl;
    QDP_abort(1);
  };

  // Get the lattice size.
  const multi1d<int>& latt_size = Layout::lattSize();
  
  // resize output array appropriately
  pion_corr_fn.resize(NUM_STAG_PIONS, latt_size[Nd-1]);

  // Correlation functions before spatial sum
  LatticeComplex corr_fn;

  // Machinery to do timeslice sums with 
  UnorderedSet timeslice;
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
  corr_fn = trace(adj(quark_props[0])*quark_props[0]);

  // Slice Sum
  pion_corr_fn[ pion_index ] = sumMulti(corr_fn, timeslice);
  pion_index++;


  // Array to describe shifts in cube
  multi1d<int> delta(Nd);

  // One link spatial pions.
  for(mu=0; mu<Nd-1; mu++) {
 
    delta = 0;
    delta[mu] = 1;
      
    corr_fn =  beta(mu)*trace(shiftDeltaProp(delta,quark_props[0])
			     *adj(quark_props[ deltaToPropIndex(delta) ]));
    
    pion_corr_fn[ pion_index ] = sumMulti(corr_fn, timeslice);
    pion_index++;
  }
    
  // One link temporal
  corr_fn = -  alpha(Nd-1)*trace(adj(quark_props[0])*quark_props[0]);
  pion_corr_fn[ pion_index ] = sumMulti(corr_fn, timeslice);
  pion_index++;

  // Two link spatial
  for(mu=0; mu<Nd-1; mu++) { 
    for(nu=mu+1; nu <Nd-1; nu++) { 
      delta = 0;
      delta[mu] = 1;
      delta[nu] = 1;

      corr_fn = - beta(mu)* beta(nu)
	                      *trace(adj(shiftDeltaProp(delta,quark_props[0]))
				     *quark_props[ deltaToPropIndex(delta) ]);
    
      pion_corr_fn[ pion_index ] = sumMulti(corr_fn, timeslice);
      pion_index++;
    }
  }

  // Two link temporal
  for(mu=0; mu<Nd-1; mu++) { 
      delta = 0;
      delta[mu] = 1;
    
      corr_fn = - beta(mu)*  alpha(Nd-1)
	                      *trace(adj(shiftDeltaProp(delta,quark_props[0]))
				     *quark_props[ deltaToPropIndex(delta) ]);
    
      pion_corr_fn[ pion_index ] = sumMulti(corr_fn, timeslice);
      pion_index++;
  }

  // Three link spatial
  for(mu=0; mu<Nd-1; mu++) { 
    for(nu=mu+1; nu <Nd-1; nu++) { 
      for(rho=nu+1; rho < Nd-1; rho++) { 

	delta = 0;
	delta[mu] = 1;
	delta[nu] = 1;
	delta[rho] = 1;

	
	corr_fn = - beta(mu) * beta(nu)* beta(rho)
	  *trace(adj(shiftDeltaProp(delta,quark_props[0]))
		 *quark_props[ deltaToPropIndex(delta) ]);
	
	pion_corr_fn[ pion_index ] = sumMulti(corr_fn, timeslice);
	pion_index++;
      }
    }
  }
  
  // Three link temporal 
  for(mu=0; mu<Nd-1; mu++) { 
    for(nu=mu+1; nu <Nd-1; nu++) { 

	delta = 0;
	delta[mu] = 1;
	delta[nu] = 1;

	corr_fn =  beta(mu)* beta(nu)*  alpha(Nd-1)
	  *trace(adj(shiftDeltaProp(delta,quark_props[0]))
		 *quark_props[ deltaToPropIndex(delta) ]);
	
	pion_corr_fn[ pion_index ] = sumMulti(corr_fn, timeslice);
	pion_index++;
    }
  }

  // Four link temporal
  delta = 0;
  delta[0] = delta[1] = delta[2] = 1;
  corr_fn = - alpha(3)* beta(0)* beta(1) * beta(2)
    *trace(adj(shiftDeltaProp(delta, quark_props[0]))
	   *quark_props[ deltaToPropIndex(delta) ] );
  
  pion_corr_fn[ pion_index ] = sumMulti(corr_fn, timeslice);
  pion_index++;

  if( pion_index != NUM_STAG_PIONS) { 
    QDPIO::cerr << "Panic! Panic! Something has gone horribly wrong" << endl;
    QDP_abort(1);
  }
}



/*! This function, converts a set of shifts stored in delta, into
  an index for the array of propagators. */
int deltaToPropIndex(multi1d<int>& delta) 
{

  if( Nd != 4 ) { 
    QDPIO::cerr << "This routine only works for Nd=4, not Nd="<<Nd << endl;
    QDP_abort(1);
  }

  /* Paranoia -- Delta has to have Nd elements */
  if( delta.size() != Nd ) { 
    QDPIO::cerr << "Delta has to have Nd elements as opposed to " << delta.size();
    QDP_abort(1);
  }

  /* Paranoia2 -- Delta can contain only 0 (no shift) and 1s (shift) */
  for( int mu = 0; mu < Nd; mu++) { 
    if( delta[mu] != 0 && delta[mu] != 1 ) { 
      QDPIO::cerr << "Delta must be made up of either zeros or ones. Element " << mu << " has value " << delta[mu] << endl;
      QDP_abort(1);
    }
  }

  /* The delta array is turned into a lexicographic index within
     the 4d hypercube */
  return(delta[0] + 2*( delta[1] + 2*(delta[2] + 2*delta[3])));
}

/*! Given an array of forward shifts (up to 1 in each dimension) and a propagator. This routine will carry out up to 4 shifts on the input, possibly 1 in each dimension */
LatticePropagator shiftDeltaProp(multi1d<int>& delta, const LatticePropagator& src)
{

  int mu;

  if( delta.size() != Nd ) { 
    QDPIO::cerr << "Delta has to have Nd elements as opposed to " << delta.size();
    QDP_abort(1);
  }

  for( mu = 0; mu < Nd; mu++) { 
    if( delta[mu] != 0 && delta[mu] != 1 ) { 
      QDPIO::cerr << "Delta must be made up of either zeros or ones. Element " << mu << " has value " << delta[mu] << endl;
      QDP_abort(1);
    }
  }

  LatticePropagator ret_val = src;
  LatticePropagator tmp;

  for( mu = 0; mu < Nd; mu++) { 
    if( delta[mu] == 1 ) { 
      // This at the moment cannot occur without a temporary
      tmp = shift(ret_val, FORWARD, mu);
      ret_val = tmp;
    }
  }

  return ret_val;
}
