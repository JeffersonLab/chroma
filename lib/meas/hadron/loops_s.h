#ifndef LOOP_S_H
#define LOOP_S_H


#include "chromabase.h"
#include "meas/hadron/stoch_var.h"
#include "meas/hadron/stag_propShift_s.h"

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

  staggered_loops(int t_len, int t_sample, 
		  multi1d<LatticeColorMatrix> & uin) : t_length(t_len) ,
						       no_sample(t_sample)
    {
      corr_fn.resize(no_sample, t_length);
      corr_fn = zero ; 
      corr.resize(t_length);
      corr = zero ; 

      u.resize(4) ; 
      if( uin.size() != 4 ) { 
	QDPIO::cerr << "staggered_hadron_corr: input guage config has wrong number of dimensions " << uin.size() << endl;
	QDP_abort(1);
      };

      u = uin ; 
      type_of_shift = GAUGE_INVAR ; 

    }

  virtual ~staggered_loops()
    {
      corr_fn.resize(1, 1);
    }



  LatticeStaggeredPropagator shift_deltaProp(multi1d<int>& delta, 
					     const 
					     LatticeStaggeredPropagator& src)

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


  void use_gauge_invar() { type_of_shift = GAUGE_INVAR ; } 
  void use_NON_gauge_invar() { type_of_shift = NON_GAUGE_INVAR ; } 


protected:

  multi2d<DComplex> corr_fn ; 
  multi1d<DComplex> corr ;
  
  string outer_tag ; 
  string inner_tag ; 
  multi1d<LatticeColorMatrix> u ; // this should handle or state

private :
  int t_length ; 
  int no_sample ; 

  Stag_shift_option type_of_shift; 

} ;




#endif
