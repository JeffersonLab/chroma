// $Id: pseudoscalar_loops_s.h,v 1.7 2005-06-27 20:22:56 mcneile Exp $
#ifndef PSEUDOSCALAR_LOOPS_S_H
#define PSEUDOSCALAR_LOOPS_S_H

#include "meas/hadron/loops_s.h"

namespace Chroma {

  class staggered_loops ; 

  class threelink_pseudoscalar_loop  : public staggered_loops
  {
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    threelink_pseudoscalar_loop(int t_len, 
				int nsample, 
				multi1d<LatticeColorMatrix> & uin,
				Stag_shift_option type_of_shift_in)  
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {
	outer_tag = "loop_gamma3gamma5_cross_one"  ; 
	inner_tag = "loop" ; 
      }

    ~threelink_pseudoscalar_loop()
      {
      }


  protected:


  } ; 



  class fourlink_pseudoscalar_loop  : public staggered_loops
  {
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    fourlink_pseudoscalar_loop(int t_len, int nsample,
			       multi1d<LatticeColorMatrix> & uin,
			       Stag_shift_option type_of_shift_in)  
      : staggered_loops(t_len,nsample,uin,type_of_shift_in)
      {

	outer_tag = "loop_gamma5_cross_one"  ; 
	inner_tag = "loop" ; 
      }

    virtual ~fourlink_pseudoscalar_loop()
      {
      }


  protected:


  } ; 

}  // end namespace Chroma

#endif
