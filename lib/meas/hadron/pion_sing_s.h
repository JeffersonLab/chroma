#ifndef PION_SING_S_H
#define PION_SING_S_H

#define NUM_STAG_PROPS   8
#define NUM_STAG_SING_PIONS   1


#include "chroma.h"

class staggered_hadron_corr ; 

class staggered_pion_singlet  : public staggered_hadron_corr
{


  public :

    void
    compute(multi1d<LatticeStaggeredPropagator>& quark_props,
		      int j_decay);


  staggered_pion_singlet(int t_len, multi1d<LatticeColorMatrix> & uin)  
   : staggered_hadron_corr(t_len,no_pion_sings,uin)
    {
      outer_tag = "SingletPseudoscalar"  ; 
      inner_tag = "Pi" ; 

      tag_names.resize(no_pion_sings) ; 

    }

  ~staggered_pion_singlet()
    {
    }


 protected:

  private :
    static const int no_pion_sings = 1 ; 


} ; 


#endif
