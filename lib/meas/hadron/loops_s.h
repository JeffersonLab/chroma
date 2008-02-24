// $Id: loops_s.h,v 3.2 2008-02-24 11:29:35 mcneile Exp $

#ifndef LOOP_S_H
#define LOOP_S_H


#include "chromabase.h"
#include "meas/hadron/stoch_var.h"
#include "meas/hadron/stag_propShift_s.h"

namespace Chroma {

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
    void dump(XMLWriter &xml_out)
    {

      multi1d<Real64> sig_sc0(t_length), imsig_sc0(t_length);

      stoch_var(corr, corr_fn, sig_sc0, imsig_sc0, 
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
      push(xml_out, "Mean");
      write(xml_out, inner_tag, corr);
      pop(xml_out);
      push(xml_out, "MeanError");
      write(xml_out, inner_tag, sig_sc0);
      pop(xml_out);
      pop(xml_out);

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
      write(xml_out, inner_tag, corr_fn[i]);
      pop(xml_out);
      pop(xml_out);
   }


 
    staggered_loops(int t_len, int t_sample, 
		    const multi1d<LatticeColorMatrix> & uin, 
		    Stag_shift_option type_of_shift_in) : 
      t_length(t_len) ,
      no_sample(t_sample) ,
      type_of_shift(type_of_shift_in) 
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

    }

    virtual ~staggered_loops()
    {
      corr_fn.resize(1, 1);
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


    // write the data in binary format 
    void binary_dump(std::string start_name)
      {

	string filename ; 
	filename = start_name + outer_tag + inner_tag ; 
	const int magic_number = 66618 ; 

	BinaryFileWriter speedy ;
	speedy.open(filename);
	write(speedy,magic_number) ;
	write(speedy,t_length) ;
	write(speedy,no_sample) ;
	write(speedy,corr_fn) ;
	write(speedy,corr) ;
	speedy.close();

      }


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


}  // end namespace Chroma


#endif
