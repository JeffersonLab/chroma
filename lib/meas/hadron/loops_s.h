#ifndef LOOP_S_H
#define LOOP_S_H


#include "chromabase.h"
#include "meas/hadron/stoch_var.h"

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
  void dump(XMLFileWriter &xml_out)
    {

      multi1d<Real64> sig_sc0(t_length), imsig_sc0(t_length);

      stoch_var(corr, corr_fn, sig_sc0, imsig_sc0, 
		t_length, no_sample);

      push(xml_out, outer_tag);
      write(xml_out, inner_tag, corr);
      pop(xml_out);




    }

  staggered_loops(int t_len, int t_sample) : t_length(t_len) ,
    no_sample(t_sample)
    {
      corr_fn.resize(no_sample, t_length);
      corr_fn = zero ; 
      corr.resize(t_length);
      corr = zero ; 
    }

  ~staggered_loops()
    {
      corr_fn.resize(1, 1);
    }


 protected:

  multi2d<DComplex> corr_fn ; 
  multi1d<DComplex> corr ;
  
  string outer_tag ; 
  string inner_tag ; 

  private :
    int no_sample ; 
  int t_length ; 



} ;




#endif
