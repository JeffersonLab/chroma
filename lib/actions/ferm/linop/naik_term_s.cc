/*  $Id: naik_term_s.cc,v 1.6 2004-12-08 21:50:24 edwards Exp $  */

/* NAIK_LINKS */

/* Construct the "triple" links */
/* used in the staggered "asqtad" action */

/* NOTE: the staggered phase factors are assumed to be included */
/*       in the gauge fields u */

/* Arguments: */

/*  u  -- gauge field (Read) */
/*  ut -- "triple-link" gauge field (Write) */

//
//  recoded by mcneile
//
//  rewritten by steve to calculate u_triple.
//

#include "actions/ferm/linop/improvement_terms_s.h"

using namespace QDP;


//
//  u(Nd), probably should check
//  ut(Nd),
//

void Triple_Links(multi1d<LatticeColorMatrix> & u,
		multi1d<LatticeColorMatrix> & ut,
		Real u0)
{
  START_CODE();
  
  LatticeColorMatrix tmp_0;
  LatticeColorMatrix tmp_1;
  Real c_3;

  int mu;

  // SZIN parameters from macros/primitives.mh
  // probably should be checked
  //
  //  enum SZINdir { BACKWARD= 0 ,  FORWARD = 1 }  ;

  
  if (Nd == 4)
  {
    c_3 = (Real)(-1) / (u0*u0*(Real)(24)); 
    QDPIO::cout << "c_3 = " << c_3 << endl;
                    
    for(mu=0; mu < Nd; ++mu)
    {
      tmp_0 = u[mu]*shift(u[mu], FORWARD, mu);
      tmp_1 = u[mu]*shift(tmp_0, FORWARD, mu);
      ut[mu] = tmp_1 * c_3;
    }
  }

  END_CODE();
}
