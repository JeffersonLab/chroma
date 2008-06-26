#ifndef PION_LOCAL_S_H
#define PION_LOCAL_S_H



#include "meas/hadron/hadron_corr_s.h"
#include "meas/hadron/stag_propShift_s.h"
#include "util/ft/sftmom.h"

namespace Chroma {

/*

This is code to compute the local pseudoscalars 
using staggered fermions.

*/

  class staggered_hadron_corr ; 

  class staggered_local_pion  : public staggered_hadron_corr
  {


  public :
    void compute(LatticeStaggeredPropagator& quark_prop_A,
		 LatticeStaggeredPropagator& quark_prop_B,
		 int j_decay) ; 

  // This is very ugly
   void compute(multi1d<LatticeStaggeredPropagator>& quark_props,
			 int j_decay) { QDP_abort(1);   }  




  void 
  compute_and_dump(
				LatticeStaggeredPropagator& quark_prop_A,
				LatticeStaggeredPropagator& quark_prop_B,
				int j_decay, int , 
				const SftMom& phases, XMLWriter& xml) ;

  
  staggered_local_pion(int t_len, 
		       const multi1d<LatticeColorMatrix> & uin,
		       Stag_shift_option type_of_shift_in = SYM_GAUGE_INVAR)
    : staggered_hadron_corr(t_len,no_pions,uin,type_of_shift_in)
    {
      outer_tag = "Pseudoscalar"  ; 
      inner_tag = "Pi" ; 
      
      tag_names.resize(no_pions) ; 
      for(int i = 0 ; i <  no_pions ; ++i ) 
	{
	  ostringstream tag;
	  tag << "re_pion" << i;
	  tag_names[i] = tag.str() ; 
	}

    }

  virtual ~staggered_local_pion()
    {
    }


  protected:

  private :
    static const int no_pions = 1 ; 

  } ;


}  // end namespace Chroma


#endif
