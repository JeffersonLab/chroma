#ifndef PSEUDOSCALAR_LOOPS_S_H
#define PSEUDOSCALAR_LOOPS_S_H


#include "meas/hadron/loops_s.h"

class staggered_loops ; 

class threelink_pseudoscalar_loop  : public staggered_loops
{
  public :


    void compute(LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi, int isample) ; 

  threelink_pseudoscalar_loop(int t_len, int nsample)  
    : staggered_loops(t_len,nsample)
    {
      outer_tag = "loop_gamma5_cross_one"  ; 
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

  fourlink_pseudoscalar_loop(int t_len, int nsample)  
    : staggered_loops(t_len,nsample)
    {
      outer_tag = "loop_gamma3gamma5_cross_one"  ; 
      inner_tag = "loop" ; 
    }

  ~fourlink_pseudoscalar_loop()
    {
    }


 protected:


} ; 


#endif
