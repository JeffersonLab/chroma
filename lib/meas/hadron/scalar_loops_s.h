#ifndef SCALAR_LOOPS_S_H
#define SCALAR_LOOPS_S_H


#include "chroma.h"

class staggered_loops ; 

class local_scalar_loop  : public staggered_loops
{


  public :


void compute(LatticeStaggeredFermion & q_source, 
	     LatticeStaggeredFermion & psi, int isample) ; 

  local_scalar_loop(int t_len, int nsample)  
    : staggered_loops(t_len,nsample)
    {
      outer_tag = "loop_one_cross_one"  ; 
      inner_tag = "loop" ; 
    }

  ~local_scalar_loop()
    {
    }


 protected:


} ; 


#endif
