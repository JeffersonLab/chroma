#ifndef LOOP_S_H
#define LOOP_S_H


#include "chroma.h"

/*

This is generic code to compute staggered bubbles.

*/


class staggered_loops
{

  public :
    virtual void compute(LatticeStaggeredFermion & q_source, 
		    LatticeStaggeredFermion & psi, int isample) = 0 ;

  /*
    Write the correlators out
  */
  void dump(int t_source, XMLFileWriter &xml_out)
    {

      multi1d<Real64> sig_sc0(t_length), imsig_sc0(t_length);
      multi1d<DComplex> sca0_loop(t_length) ; 
#ifdef NNNNNN
      stoch_var(sca0_loop, corr_fn, sig_sc0, imsig_sc0, 
		t_length, Nsamp);
#endif

      push(xml_out, tag_name);
      write(xml_out, inner_tag, sca0_loop);
      pop(xml_out);




    }

  staggered_loops(int t_len, int t_sample) : t_length(t_len) ,
    no_sample(t_sample)
    {
      corr_fn.resize(no_sample, t_length);
    }

  ~staggered_loops()
    {
      corr_fn.resize(1, 1);
    }


 protected:

  multi2d<DComplex> corr_fn ; 
  string outer_tag ; 
  string inner_tag ; 
  string tag_name ; 

  private :
    int no_sample ; 
  int t_length ; 



} ;




#endif
