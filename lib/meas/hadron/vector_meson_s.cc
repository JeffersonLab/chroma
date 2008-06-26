/*! File:  vector_meson_s.cc   
 *
 * The routines in this file computes some of the staggered rhos
 *
 * Only the local vectors are computed.

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

/*! StaggeredVectorMesons
 *
 * This routine computes the three local vector mesons.
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


  // Also for now we want j_decay to be 3.
  switch( j_decay ) { 
  case 3:
    break;
    
  default:
    QDPIO::cerr << "staggeredVectors: j_decay must be 3 for just now. It is " << j_decay << endl;
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
  
  // Counters/Indices
  int sca_index = 0;

  //
  //
  //
    corr_fn_s =  StagPhases::alpha(1)*
        trace(quark_props[0]*adj(quark_props[ 0 ]));
    corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
    tag_names[sca_index] = "gamma_x_CROSS_gamma_x" ; 
    sca_index++;
      
    corr_fn_s =  StagPhases::alpha(1)*StagPhases::alpha(2)*
                trace(quark_props[0]*adj(quark_props[0]));

    corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
    tag_names[sca_index] = "gamma_y_CROSS_gamma_y" ; 
    sca_index++;


    corr_fn_s =  StagPhases::alpha(2)*StagPhases::alpha(3)*
                             trace(quark_props[0]*adj(quark_props[0] ));
    corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
    tag_names[sca_index] = "gamma_z_CROSS_gamma_z" ; 
    sca_index++;


    tag_names[sca_index] = "M_VT" ; 

    for(int t= 0 ;  t < latt_size[Nd-1] ; ++t)
    {
      corr_fn[sca_index][t] = ( corr_fn[0][t] 
				+ corr_fn[1][t] 
				+ corr_fn[2][t] ) / 3.0 ; 
    }


}


/**
     Compute the local vector meson correlators using a single
     local quark propagator.

 * I have checked this routine against the MILC code. 

**/


void 
vector_meson::compute(
   LatticeStaggeredPropagator & quark_props,
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
    QDPIO::cerr << "staggeredVectors: j_decay must be 3 for just now. It is " << j_decay << endl;
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
  
  // Counters/Indices
  int sca_index = 0;

  //
  //
  //

  corr_fn_s =  StagPhases::alpha(1)*
    trace(quark_props*adj(quark_props));
  corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
  tag_names[sca_index] = "gamma_x_CROSS_gamma_x" ; 
  sca_index++;

    corr_fn_s =  StagPhases::alpha(1)*StagPhases::alpha(2)*
                trace(quark_props*adj(quark_props));

    corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
    tag_names[sca_index] = "gamma_y_CROSS_gamma_y" ; 
    sca_index++;


    corr_fn_s =  StagPhases::alpha(2)*StagPhases::alpha(3)*
                             trace(quark_props*adj(quark_props ));
    corr_fn[ sca_index ] = sumMulti(corr_fn_s, timeslice);
    tag_names[sca_index] = "gamma_z_CROSS_gamma_z" ; 
    sca_index++;


    tag_names[sca_index] = "M_VT" ; 

    for(int t= 0 ;  t < latt_size[Nd-1] ; ++t)
    {
      corr_fn[sca_index][t] = (   corr_fn[0][t] 
				+ corr_fn[1][t] 
				+ corr_fn[2][t] ) / 3.0 ; 
    }


}

}  // end namespace Chroma
