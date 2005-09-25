// $Id: pion_sing_s.h,v 2.0 2005-09-25 21:04:35 edwards Exp $
#ifndef PION_SING_S_H
#define PION_SING_S_H

#include "chromabase.h"
#include "meas/hadron/hadron_corr_s.h"

namespace Chroma {


  class staggered_hadron_corr ; 

  class staggered_pion_singlet  : public staggered_hadron_corr
  {


  public :

    void compute(
      LatticeStaggeredPropagator local_quark_prop,
      LatticeStaggeredPropagator four_shift_quark_prop,
      int j_decay) ;


    void
    compute(multi1d<LatticeStaggeredPropagator>& quark_props,
	    int j_decay) { } 


    staggered_pion_singlet(int t_len, 
			   const multi1d<LatticeColorMatrix> & uin,
			   Stag_shift_option type_of_shift_in = SYM_GAUGE_INVAR)  
      : staggered_hadron_corr(t_len,no_pion_sings,uin,type_of_shift_in)
      {
	outer_tag = "SingletPseudoscalar"  ; 
	inner_tag = "Pi" ; 

	tag_names.resize(no_pion_sings) ; 

      }

    virtual ~staggered_pion_singlet()
      {
      }


  protected:

  private :
    static const int no_pion_sings = 1 ; 


  } ; 

}  // end namespace Chroma

#endif
