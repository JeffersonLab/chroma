/*! File:  vector_meson_s.cc   
 *
 * The routines in this file computes some of the staggered rhos
 * (This routine has not yet been CHECKED !!!!!)
 * 
 * BEWARE: These routines ASSUME that the quark propagators have been
 * calculated in the Coulomb (or other spatially fixed) gauge.
 * It is not gauge invariant. It could be made to be so by adding some
 * parallel transport, however folklore claims that increases noise
 *
 * YOU HAVE BEEN WARNED.
 */

#include "meas/hadron/vector_meson_s.h"
#include "meas/hadron/stag_propShift_s.h"
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
 *
 *       j_decay          -- The time direction (has to be Nd-1 for now)
 *                           (Read)
 */
 

// THis brings the staggered phases alpha and beta into the namespace

void 
vector_meson::compute(
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

  // 
  // Array to describe shifts in cube
  multi1d<int> delta(Nd);


  //
  //
  //

    mu = 0 ; 
    delta = 0;
    delta[mu] = 1;
      
    corr_fn_s =  StagPhases::alpha(1)*trace(shift_deltaProp(delta,quark_props[0])
                             *adj(quark_props[ deltaToPropIndex(delta) ]));
    corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
    tag_names[sca_index] = "gamma_x_CROSS_gamma_x" ; 
    sca_index++;


    mu = 1 ; 
    delta = 0;
    delta[mu] = 1;
      
    corr_fn_s =  StagPhases::alpha(1)*trace(shift_deltaProp(delta,quark_props[0])
                             *adj(quark_props[ deltaToPropIndex(delta) ]));
    corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
    tag_names[sca_index] = "gamma_y_CROSS_gamma_y" ; 
    sca_index++;


    mu = 2 ; 
    delta = 0;
    delta[mu] = 1;
      
    corr_fn_s =  StagPhases::alpha(1)*trace(shift_deltaProp(delta,quark_props[0])
                             *adj(quark_props[ deltaToPropIndex(delta) ]));
    corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
    tag_names[sca_index] = "gamma_z_CROSS_gamma_z" ; 
    sca_index++;

    //
    //
    //

#ifdef NNNNNNNNNNNNNNn
    corr_fn[ sca_index ] = ( 
			    corr_fn[ sca_index ] 
			    + corr_fn[ sca_index ] 
			    + corr_fn[ sca_index ] ) / 3.0 ; 
#endif

    tag_names[sca_index] = "M_VT" ; 


}

}  // end namespace Chroma
