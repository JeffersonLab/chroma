#ifndef PIONS_S_H
#define PIONS_S_H

#define NUM_STAG_PROPS   8
#define NUM_STAG_PIONS   16


#include "chroma.h"

/*

 This is code to compute all 16 pseudoscalars 
 using staggered fermions.

*/

class staggered_hadron_corr ; 

class staggered_pions  : public staggered_hadron_corr
{


  public :
    void compute(multi1d<LatticeStaggeredPropagator>& quark_props,
	  int j_decay) ;


  staggered_pions(int t_len)  : staggered_hadron_corr(t_len,no_pions)
    {
      outer_tag = "Here_are_all_16_pions"  ; 
      inner_tag = "Pi" ; 

      tag_names.resize(no_pions) ; 
      for(int i = 0 ; i <  no_pions ; ++i ) 
	{
	  ostringstream tag;
	  tag << "re_pion" << i;
	  tag_names[i] = tag.str() ; 
	}

    }

  ~staggered_pions()
    {
    }


 protected:

  private :
    static const int no_pions = 16 ; 

} ;




#endif
