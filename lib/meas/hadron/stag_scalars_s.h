#ifndef STAG_SCALARS_S_H
#define STAG_SCALARS_S_H

#define NUM_STAG_PROPS   8
#define NUM_STAG_PIONS   16


#include "meas/hadron/hadron_corr_s.h"

namespace Chroma {

  class staggered_hadron_corr ; 

  class staggered_scalars  : public staggered_hadron_corr
  {


  public :

    void
    compute(multi1d<LatticeStaggeredPropagator>& quark_props,
	    int j_decay);


    staggered_scalars(int t_len, 
		      const multi1d<LatticeColorMatrix> & uin,
		      Stag_shift_option type_of_shift_in = SYM_GAUGE_INVAR)  
      : staggered_hadron_corr(t_len,no_scalar,uin,type_of_shift_in)
      {
	outer_tag = "Here_are_all_16_scalars"  ; 
	inner_tag = "Sc" ; 

	tag_names.resize(no_scalar) ; 
	for(int i = 0 ; i <  no_scalar ; ++i ) 
	{
	  ostringstream tag;
	  tag << "re_sc" << i;
	  tag_names[i] = tag.str() ; 
	}

      }

    virtual ~staggered_scalars()
      {
      }


  protected:

  private :
    static const int no_scalar = 16 ; 


  } ; 

}  // end namespace Chroma

#endif
