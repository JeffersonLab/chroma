#ifndef HADRON_CORR_S_H
#define HADRON_CORR_S_H

#define NUM_STAG_PROPS   8
#define NUM_STAG_PIONS   16


#include "chroma.h"

/*

This is generic code to compute staggered correlators.

*/


class staggered_hadron_corr
{

  public :
    virtual void compute(multi1d<LatticeStaggeredPropagator>& quark_props,
	  int j_decay) = 0 ;

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
      push(xml_out, outer_tag);
      for(int i=0; i < no_channel; i++) {
	Pi = re_corr_fn[i];
	//	ostringstream tag;
	//tag << "re_pion" << i;
	// push(xml_out, tag.str());
	push(xml_out, tag_names[i]);
	write(xml_out, inner_tag, Pi);
	pop(xml_out);
      }
      pop(xml_out);



    }

  staggered_hadron_corr(int t_len, int t_chan) : t_length(t_len) ,
    no_channel(t_chan)
    {
      corr_fn.resize(no_channel, t_length);
    }

  ~staggered_hadron_corr()
    {
      corr_fn.resize(1, 1);
    }


 protected:

  multi2d<DComplex> corr_fn ; 
  string outer_tag ; 
  string inner_tag ; 
  multi1d<string> tag_names ; 

  private :
    int no_channel ; 
  int t_length ; 



} ;




#endif
