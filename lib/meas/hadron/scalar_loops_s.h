// $Id: scalar_loops_s.h,v 1.6 2005-01-14 18:42:36 edwards Exp $
#ifndef SCALAR_LOOPS_S_H
#define SCALAR_LOOPS_S_H

#include "meas/hadron/loops_s.h"

namespace Chroma {

  class staggered_loops ; 

  class local_scalar_loop  : public staggered_loops
  {
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    local_scalar_loop(int t_len, int nsample,
		      multi1d<LatticeColorMatrix> & uin)  
      : staggered_loops(t_len,nsample,uin)
      {
	outer_tag = "loop_one_cross_one"  ; 
	inner_tag = "loop" ; 
      }

    virtual ~local_scalar_loop()
      {
      }


  protected:


  } ; 



  class non_local_scalar_loop  : public staggered_loops
  {
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

    non_local_scalar_loop(int t_len, int nsample,
			  multi1d<LatticeColorMatrix> & uin)  
      : staggered_loops(t_len,nsample,uin)
      {
	outer_tag = "loop_gamma3_cross_one"  ; 
	inner_tag = "loop" ; 
      }

    virtual ~non_local_scalar_loop()
      {
      }


  protected:


  } ; 

}  // end namespace Chroma

#endif
