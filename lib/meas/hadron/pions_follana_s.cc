
#include "chroma.h"
#include "mesphas_follana_s.h"
#include "pions_follana_s.h"

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


int deltaToPropIndex(multi1d<int>& delta) 
{
  if( delta.size() != Nd ) { 
    QDPIO::cerr << "Delta has to have Nd elements as opposed to " << delta.size();
    QDP_abort(1);
  }

  for( int mu = 0; mu < Nd; mu++) { 
    if( delta[mu] != 0 && delta[mu] != 1 ) { 
      QDPIO::cerr << "Delta must be made up of either zeros or ones. Element " << mu << " has value " << delta[mu] << endl;
      QDP_abort(1);
    }
  }

  return(delta[0] + 2*( delta[1] + 2*(delta[2] + 2*delta[3])));
}

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
      tmp = shift(ret_val, FORWARD, mu);
      ret_val = tmp;
    }
  }

  return ret_val;
}

 



void 
staggeredPionsFollana(multi1d<LatticePropagator>& quark_props,
		      multi2d<DComplex>& pion_corr_fn,
		      int j_decay)
{
  if( Nd != 4 ) { 
    QDPIO::cerr << "The no of dimensions should be 4 for now. It is: " << Nd << endl;
    QDP_abort(1);
  }

  // Some dumb things to check -- that quark_props has 8
  // members
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
  
  // Make the array the right size
  pion_corr_fn.resize(NUM_STAG_PIONS, latt_size[Nd-1]);

  // Correlation functions before spatial sum
  LatticeComplex corr_fn;

  // Machinery to do timeslice sums with 
  UnorderedSet timeslice;
  timeslice.make(TimeSliceFunc(Nd-1));

  multi1d<LatticeInteger> alpha(Nd);
  multi1d<LatticeInteger> beta(Nd);

  // Get the phases
  mesPhasFollana(alpha, beta);
  
  int pion_index = 0;
  int i;
  int mu;
  int nu;
  int rho;

  // Goldstone Pion
  corr_fn = trace(adj(quark_props[0])*quark_props[0]);
  pion_corr_fn[ pion_index ] = sumMulti(corr_fn, timeslice);
  pion_index++;

  multi1d<int> delta(Nd);

  // One link spatial
  for(mu=0; mu<Nd-1; mu++) { 
    delta = 0;
    delta[mu] = 1;
      
    corr_fn = beta[mu]*trace(shiftDeltaProp(delta,quark_props[0])
			     *adj(quark_props[ deltaToPropIndex(delta) ]));
    
    pion_corr_fn[ pion_index ] = sumMulti(corr_fn, timeslice);
    pion_index++;
  }
    
  // One link temporal
  corr_fn = -alpha[Nd-1]*trace(adj(quark_props[0])*quark_props[0]);
  pion_corr_fn[ pion_index ] = sumMulti(corr_fn, timeslice);
  pion_index++;

  // Two link spatial
  for(mu=0; mu<Nd-1; mu++) { 
    for(nu=mu+1; nu <Nd-1; nu++) { 
      delta = 0;
      delta[mu] = 1;
      delta[nu] = 1;

      corr_fn = -beta[mu]*beta[nu]
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
    
      corr_fn = -beta[mu]*alpha[Nd-1]
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

	
	corr_fn = -beta[mu]*beta[nu]*beta[rho]
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

	corr_fn = beta[mu]*beta[nu]*alpha[Nd-1]
	  *trace(adj(shiftDeltaProp(delta,quark_props[0]))
		 *quark_props[ deltaToPropIndex(delta) ]);
	
	pion_corr_fn[ pion_index ] = sumMulti(corr_fn, timeslice);
	pion_index++;
    }
  }

  // Four link temporal
  delta = 0;
  delta[0] = delta[1] = delta[2] = 1;
  corr_fn = -alpha[3]*beta[0]*beta[1]*beta[2]
    *trace(adj(shiftDeltaProp(delta, quark_props[0]))
	   *quark_props[ deltaToPropIndex(delta) ] );
  
  pion_corr_fn[ pion_index ] = sumMulti(corr_fn, timeslice);
  pion_index++;

  if( pion_index != NUM_STAG_PIONS) { 
    QDPIO::cerr << "Panic! Panic! Something has gone horribly wrong" << endl;
    QDP_abort(1);
  }
}
