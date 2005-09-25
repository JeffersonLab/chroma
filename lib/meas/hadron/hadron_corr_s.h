// $Id: hadron_corr_s.h,v 2.0 2005-09-25 21:04:35 edwards Exp $

#ifndef HADRON_CORR_S_H
#define HADRON_CORR_S_H

#define NUM_STAG_PROPS   8
#define NUM_STAG_PIONS   16


#include "chromabase.h"
#include "stag_propShift_s.h"

namespace Chroma 
{

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
    void dump(int t_source, XMLWriter &xml_out)
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

    staggered_hadron_corr(int t_len, int t_chan, 
			  const multi1d<LatticeColorMatrix> & uin, 
			  Stag_shift_option type_of_shift_in  ) 
      : t_length(t_len) ,
	no_channel(t_chan)
    {
      corr_fn.resize(no_channel, t_length);
      u.resize(4) ; 

      if( uin.size() != 4 ) { 
	QDPIO::cerr << "staggered_hadron_corr: input guage config has wrong number of dimensions " << uin.size() << endl;
	QDP_abort(1);
      };

      u = uin ; 
      type_of_shift = type_of_shift_in ; 
    }


    virtual ~staggered_hadron_corr()
    {
      corr_fn.resize(1, 1);
    }


  protected:

    multi2d<DComplex> corr_fn ; 
    string outer_tag ; 
    string inner_tag ; 
    multi1d<string> tag_names ; 
    multi1d<LatticeColorMatrix> u ; // this should handle or state

    LatticeStaggeredPropagator shift_deltaProp(multi1d<int>& delta, 
					       const LatticeStaggeredPropagator& src)

    {

      switch (type_of_shift)
	{
	case NON_GAUGE_INVAR :
	  return shiftDeltaProp(delta,src) ;
	  break ;
	case GAUGE_INVAR :
	  return shiftDeltaPropCov(delta,src,u,false) ;
	  break ;
	case SYM_GAUGE_INVAR :
	  return shiftDeltaPropCov(delta,src,u,true) ; // symm shifting
	  break ;
	case SYM_NON_GAUGE_INVAR:
	  return shiftDeltaProp(delta,src,true) ; // symm shifting
	  break ;
	default :
	  /**************************************************************************/
 
	  QDPIO::cerr << "Shift type " << type_of_shift << " unsupported." << endl;
	  QDP_abort(1);
	}

    }


  private :
    int t_length ; 
    int no_channel ; 

    Stag_shift_option type_of_shift; 




  } ;


}  // end namespace Chroma


#endif
