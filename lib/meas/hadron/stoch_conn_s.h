// $Id: stoch_conn_s.h,v 3.0 2006-04-03 04:59:00 edwards Exp $

#ifndef STOCH_CONN_S_H
#define STOCH_CONN_S_H


#include "chromabase.h"
#include "meas/hadron/stoch_var.h"
#include "meas/hadron/stag_propShift_s.h"

namespace Chroma {

  /*

  This is generic code to compute staggered bubbles.

  */


  class stoch_conn_corr
  {

  public :
    virtual void compute(LatticeStaggeredFermion & q_source1,
		 LatticeStaggeredFermion & q_source2,
		 LatticeStaggeredFermion & psi1,
		 LatticeStaggeredFermion & psi2,
		 int isample) = 0 ;

    /*
      Write the correlators out
    */
    void dump(XMLWriter &xml_out)
    {

      multi1d<Real64> sig_sc1(t_length), imsig_sc1(t_length);
      multi1d<Real64> sig_sc2(t_length), imsig_sc2(t_length);


      stoch_var(corr1, corr_fn1, sig_sc1, imsig_sc1, 
		t_length, no_sample);
      stoch_var(corr2, corr_fn2, sig_sc2, imsig_sc2, 
		t_length, no_sample);

      push(xml_out, outer_tag);
      switch (type_of_shift)
	{
	case NON_GAUGE_INVAR :
          write(xml_out, "SHIFT", "NON_GAUGE_INVAR");
	  break ;
	case GAUGE_INVAR :
          write(xml_out, "SHIFT", "GAUGE_INVAR");
	  break ;
	case SYM_GAUGE_INVAR :
          write(xml_out, "SHIFT", "SYM_GAUGE_INVAR");
	  break ;
	case SYM_NON_GAUGE_INVAR:
          write(xml_out, "SHIFT", "SYM_NON_GAUGE_INVAR");
	  break ;
	}

      push(xml_out, "Mean1");
      write(xml_out, inner_tag1, corr1);
      pop(xml_out);    // Mean1

      push(xml_out, "MeanError1");
      write(xml_out, inner_tag1, sig_sc1);
      pop(xml_out);    // MeanError1

      //      pop(xml_out);

      push(xml_out, "Mean2");
      write(xml_out, inner_tag2, corr2);
      pop(xml_out);     // Mean2

      push(xml_out, "MeanError2");
      write(xml_out, inner_tag2, sig_sc2);
      pop(xml_out);    //MeanError2

      pop(xml_out);   // outer_tag

    }


   void dump(XMLWriter &xml_out, int &i){

      string   tag;
      char *cnum;
      string snum;

      cnum=(char*)malloc(10*sizeof(char));

      sprintf(cnum,"%d",i);

      snum=cnum;
      free(cnum);

      string strzeros(6-snum.size(),'0');

      tag = "Meas"+strzeros+snum;

      push(xml_out, outer_tag);
      push(xml_out, tag);
      write(xml_out, inner_tag1, corr_fn1[i]);
      //      pop(xml_out);

      //      push(xml_out, tag);
      write(xml_out, inner_tag2, corr_fn2[i]);
      pop(xml_out);
      pop(xml_out);
    }

    stoch_conn_corr(int t_len, int t_sample, 
		    const multi1d<LatticeColorMatrix> & uin, 
		    Stag_shift_option type_of_shift_in) : 
      t_length(t_len) ,
      no_sample(t_sample) ,
      type_of_shift(type_of_shift_in) 
    {
      corr_fn1.resize(no_sample, t_length);
      corr_fn2.resize(no_sample, t_length);
      corr_fn1 = zero ; 
      corr_fn2 = zero ; 
      corr1.resize(t_length);
      corr2.resize(t_length);
      corr1 = zero ; 
      corr2 = zero ; 

      u.resize(4) ; 
      if( uin.size() != 4 ) { 
	QDPIO::cerr << "staggered_hadron_corr: input guage config has wrong number of dimensions " << uin.size() << endl;
	QDP_abort(1);
      };

      u = uin ; 

    }

    virtual ~stoch_conn_corr()
    {
      corr_fn1.resize(1, 1);
      corr_fn2.resize(1, 1);
    }



    LatticeStaggeredFermion shift_deltaProp(multi1d<int>& delta, 
					       const 
					       LatticeStaggeredFermion& src)

    {

      switch (type_of_shift)
	{
	case NON_GAUGE_INVAR :
          cout << "ERROR SHIFT: GI\n"  ; exit(0) ;
	  //	  return shiftDeltaProp(delta,src) ;
	  break ;
	case GAUGE_INVAR :
          cout << "ERROR SHIFT: GI\n"  ; exit(0) ;
	  //	  return shiftDeltaPropCov(delta,src,u,false) ;
	  break ;
	case SYM_GAUGE_INVAR :
	  return shiftDeltaPropCov(delta,src,u,true) ; // symm shifting
	  break ;
	case SYM_NON_GAUGE_INVAR:
          cout << "ERROR SHIFT: SNGI\n"  ; exit(0) ;
	  //	  return shiftDeltaProp(delta,src,true) ; // symm shifting
	  break ;
	default :
	  /**************************************************************************/
 
	  QDPIO::cerr << "Shift type " << type_of_shift << " unsupported." << endl;
	  QDP_abort(1);
	}

       return zero;   // make compiler happy
    }



  protected:

    multi2d<DComplex> corr_fn1 ; 
    multi2d<DComplex> corr_fn2 ; 
    multi1d<DComplex> corr1 ;
    multi1d<DComplex> corr2 ;
  
    string outer_tag ; 
    string inner_tag1 ; 
    string inner_tag2 ; 
    multi1d<LatticeColorMatrix> u ; // this should handle or state

  private :
    int t_length ; 
    int no_sample ; 

    Stag_shift_option type_of_shift; 

  } ;


}  // end namespace Chroma


#endif
