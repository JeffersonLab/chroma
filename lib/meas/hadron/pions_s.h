#ifndef PIONS_S_H
#define PIONS_S_H

#define NUM_STAG_PROPS   8
#define NUM_STAG_PIONS   16


#include "chroma.h"

/*

 This is code to compute all 16 pseudoscalars 
 using staggered fermions.

*/


class staggered_pions
{


  public :
    void compute(multi1d<LatticeStaggeredPropagator>& quark_props,
	  int j_decay) ;

  /*
    Write the correlators out
  */
  void dump(int t_source, XMLFileWriter &xml_out)
    {
      multi2d<Real> pi_corr(no_channel, t_length);
      multi2d<Real> re_corr_fn(no_channel, t_length);
      pi_corr = zero;

      // Take the real part of the correlator and average over the sources
      for(int i=0; i < no_channel ; i++){
	for(int t=0; t < t_length; t++){
	  int t_eff = (t - t_source + t_length)% t_length;
	  re_corr_fn[i][t_eff] = real(corr_fn[i][t]);
	}
      }
                  
      multi1d<Real> Pi(t_length) ; 
      push(xml_out, "Here_are_all_16_pions");
      for(int i=0; i < no_channel; i++) {
	Pi = re_corr_fn[i];
	ostringstream tag;
	tag << "re_pion" << i;
	push(xml_out, tag.str());
	write(xml_out, "Pi", Pi);
	pop(xml_out);
      }
      pop(xml_out);



    }

  staggered_pions(int t_len) 
    {
      t_length = t_len ; 
      corr_fn.resize(no_channel, t_length);
    }

  ~staggered_pions()
    {
      corr_fn.resize(1, 1);
    }


 protected:

  multi2d<DComplex> corr_fn ; 

  private :
    static const int no_channel = 16 ; 
  int t_length ; 

} ;




#endif
