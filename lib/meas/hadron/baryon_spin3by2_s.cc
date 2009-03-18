/*
   Code for j=3/2 Staggered Baryons.

Lattice Baryons With Staggered Fermions.
Maarten F.L. Golterman, Jan Smit (Amsterdam U.)
Published in Nucl.Phys.B255:328,1985

The Quenched spectrum with staggered fermions.
Rajan Gupta (Los Alamos) , Gerald Guralnik (Brown U.) , Gregory W. Kilcup (Ohio State U.) , Stephen R. Sharpe 
Published in Phys.Rev.D43:2003-2026,1991

The class 7 operators has been compared against the MILC
code. The class 4 operator is still experimental!

*/

#include "chromabase.h"

namespace Chroma {

  // some prototypes

  void baryon_s(
		LatticeStaggeredPropagator & quark_propagator_in_a,
		LatticeStaggeredPropagator & quark_propagator_in_b,
		LatticeStaggeredPropagator & quark_propagator_in_c,
		multi1d<Complex> & barprop,
		multi1d<int> & t_source,
		int j_decay, int bc_spec) ;


  //
  //  This is the quasi-local operator that is used by 
  //  MILC from the above Gupta paper
  //
  //  This is known as a class 7 operator by Smit and Golterman.
  //
  //

  void baryon_class7_s(
		LatticeStaggeredPropagator & quark_propagator_in_a,
		LatticeStaggeredPropagator & quark_propagator_in_b,
		LatticeStaggeredPropagator & quark_propagator_in_c,
		multi1d<Complex> & barprop,
		multi1d<int> & t_source,
		int j_decay, int bc_spec) 
  {
    StopWatch swatch;
    swatch.start();

    // shift the quark propagators
    LatticeStaggeredPropagator quark_propagator_in_a_shift ;
    LatticeStaggeredPropagator quark_propagator_in_b_shift ;
    LatticeStaggeredPropagator quark_propagator_in_c_shift ;

    // shift in x direction
    quark_propagator_in_a_shift = shift(quark_propagator_in_a,FORWARD,0) + 
      shift(quark_propagator_in_a,BACKWARD,0) ;

    quark_propagator_in_b_shift = shift(quark_propagator_in_b,FORWARD,1) + 
      shift(quark_propagator_in_b,BACKWARD,1) ;

    quark_propagator_in_c_shift = shift(quark_propagator_in_c,FORWARD,2) + 
      shift(quark_propagator_in_c,BACKWARD,2) ;

    // call the Wick contarction routine
    baryon_s(quark_propagator_in_a_shift,
	     quark_propagator_in_b_shift,
	     quark_propagator_in_c_shift, 
	     barprop, t_source, j_decay, bc_spec) ;  

    // 8 = 2**3 for average in the derivative 
    // 6 is the factor in the definition  of the operator

    //    barprop /= 6.0 * 8.0 ;

    // chose normalisation to agree with MILC code
    barprop /= 6.0  ;

    swatch.stop();
    double time_in_sec  = swatch.getTimeInSeconds();

    QDPIO::cout << "baryon_class7_s took " << time_in_sec
                << " secs" << endl;

  }



  //
  //  This is the quasi-local NLT operator 
  //  from the Gupta paper. One of the quark fields are
  //  shifted forward by one time slice.
  //
  //   This is known as a class 7 operator by Smit and Golterman.
  //
  //

  void baryon_class7_NLT_s(
		LatticeStaggeredPropagator & quark_propagator_in_a,
		LatticeStaggeredPropagator & quark_propagator_in_b,
		LatticeStaggeredPropagator & quark_propagator_in_c,
		multi1d<LatticeColorMatrix> & u , 
		multi1d<Complex> & barprop,
		multi1d<int> & t_source,
		int j_decay, int bc_spec) 
  {
    StopWatch swatch;
    swatch.start();

    LatticeStaggeredPropagator quark_propagator_in_a_TIME_shift ;
    // use a covariant shift in time because we only fix to Coulomb
    // gauge. This may make the correlators more noisy

    quark_propagator_in_a_TIME_shift = u[j_decay] * shift(quark_propagator_in_a_TIME_shift,FORWARD,j_decay) ; 


    // shift the quark propagators
    LatticeStaggeredPropagator quark_propagator_in_a_shift ;
    LatticeStaggeredPropagator quark_propagator_in_b_shift ;
    LatticeStaggeredPropagator quark_propagator_in_c_shift ;

    // shift in x direction
    quark_propagator_in_a_shift = 
      shift(quark_propagator_in_a_TIME_shift,FORWARD,0) + 
      shift(quark_propagator_in_a_TIME_shift,BACKWARD,0) ;

    quark_propagator_in_b_shift = shift(quark_propagator_in_b,FORWARD,1) + 
      shift(quark_propagator_in_b,BACKWARD,1) ;

    quark_propagator_in_c_shift = shift(quark_propagator_in_c,FORWARD,2) + 
      shift(quark_propagator_in_c,BACKWARD,2) ;

    // call the Wick contarction routine
    baryon_s(quark_propagator_in_a_shift,
	     quark_propagator_in_b_shift,
	     quark_propagator_in_c_shift, 
	     barprop, t_source, j_decay, bc_spec) ;  

    // 8 = 2**3 for average in the derivative 
    // 6 is the factor in the definition  of the operator

    //    barprop /= 6.0 * 8.0 ;

    // chose normalisation to agree with MILC code
    barprop /= 6.0  ;

    swatch.stop();
    double time_in_sec  = swatch.getTimeInSeconds();

    QDPIO::cout << "baryon_class7_s took " << time_in_sec
                << " secs" << endl;

  }


  /**
      This is the class 4 staggered baryon operator
      that couples to the j=3/2 state.

  **/


  void baryon_class4_s(
		LatticeStaggeredPropagator & quark_propagator_in_a,
		LatticeStaggeredPropagator & quark_propagator_in_b,
		LatticeStaggeredPropagator & quark_propagator_in_c,
		multi1d<Complex> & barprop,
		multi1d<int> & t_source,
		int j_decay, int bc_spec) 
  {
    StopWatch swatch;
    swatch.start();

    int t_length = barprop.size() ;
    multi1d<Complex>  barprop_tmp(t_length) ; 

    // shift the quark propagators
    LatticeStaggeredPropagator quark_propagator_in_a_shift ;
    LatticeStaggeredPropagator quark_propagator_in_b_shift ;
    LatticeStaggeredPropagator quark_propagator_in_c_shift ;

    LatticeStaggeredPropagator tmp ;

    // ----- shift in x direction -----
    quark_propagator_in_a_shift = shift(quark_propagator_in_a,FORWARD,0) ;

    tmp = shift(quark_propagator_in_b,FORWARD,0) ;
    quark_propagator_in_b_shift = shift(tmp,FORWARD,1) + 
      shift(tmp,BACKWARD,1) ;

    tmp = shift(quark_propagator_in_c,FORWARD,0) ;
    quark_propagator_in_c_shift = shift(tmp,FORWARD,2) + 
      shift(tmp,BACKWARD,2) ;

    // call the Wick contarction routine
    baryon_s(quark_propagator_in_a_shift,
	     quark_propagator_in_b_shift,
	     quark_propagator_in_c_shift, 
	     barprop_tmp, t_source, j_decay, bc_spec) ;  

    barprop_tmp /= (-6.0 * 4.0) ;
    barprop = barprop_tmp ;
    

    // ----- shift in y direction -----
    quark_propagator_in_a_shift = shift(quark_propagator_in_a,FORWARD,1) ;

    tmp = shift(quark_propagator_in_b,FORWARD,1) ;
    quark_propagator_in_b_shift = shift(tmp,FORWARD,0) + 
      shift(tmp,BACKWARD,0) ;

    tmp = shift(quark_propagator_in_c,FORWARD,1) ;
    quark_propagator_in_c_shift = shift(tmp,FORWARD,2) + 
      shift(tmp,BACKWARD,2) ;

    // call the Wick contarction routine
    baryon_s(quark_propagator_in_a_shift,
	     quark_propagator_in_b_shift,
	     quark_propagator_in_c_shift, 
	     barprop_tmp, t_source, j_decay, bc_spec) ;  

    barprop_tmp /= (6.0 * 4.0) ;
    barprop += barprop_tmp ;


    // ----- shift in z direction -----
    quark_propagator_in_a_shift = shift(quark_propagator_in_a,FORWARD,2) ;

    tmp = shift(quark_propagator_in_b,FORWARD,2) ;
    quark_propagator_in_b_shift = shift(tmp,FORWARD,0) + 
      shift(tmp,BACKWARD,0) ;

    tmp = shift(quark_propagator_in_c,FORWARD,2) ;
    quark_propagator_in_c_shift = shift(tmp,FORWARD,1) + 
      shift(tmp,BACKWARD,1) ;

    // call the Wick contarction routine
    baryon_s(quark_propagator_in_a_shift,
	     quark_propagator_in_b_shift,
	     quark_propagator_in_c_shift, 
	     barprop_tmp, t_source, j_decay, bc_spec) ;  

    barprop_tmp /= (-6.0 * 4.0) ;
    barprop += barprop_tmp ;

    swatch.stop();
    double time_in_sec  = swatch.getTimeInSeconds();



    QDPIO::cout << "baryon_class4_s took " << time_in_sec
                << " secs" << endl;



  }



} // end namespace Chroma

