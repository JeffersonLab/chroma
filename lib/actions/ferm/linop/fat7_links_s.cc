/*  $Id: fat7_links_s.cc,v 1.4 2004-12-12 21:22:15 edwards Exp $  */

/* FAT7_LINKS */

/* Construct the "fat" links with staples up to 7 links long */
/* used in the staggered "asqtad" action */

/* NOTE: the staggered phase factors are assumed to be included */
/*       in the gauge fields u */

/* Arguments: */

/*  u  -- gauge field (Read) */
/*  uf -- "fat-link" gauge field (Write) */

//
//  recoded by mcneile
//
//  NO CHECKERBOARDING HERE

#include "actions/ferm/linop/improvement_terms_s.h"

using namespace QDP;

namespace Chroma 
{ 

//
//  u(Nd), probably should check
//  uf(Nd),
//

void Fat7_Links(multi1d<LatticeColorMatrix> & u,
		multi1d<LatticeColorMatrix> & uf,
		Real u0)
{
  START_CODE();
  
  LatticeColorMatrix tmp_0;
  LatticeColorMatrix tmp_1;
  LatticeColorMatrix tmp_2;
  LatticeColorMatrix tmp_3;
  LatticeColorMatrix tmp_4;
  Real c_1l;
  Real c_3l;
  Real c_5l;
  Real c_7l;
  Real c_Lepage;
  int mu;
  int nu;
  int rho;
  int sigma;

  // SZIN parameters from macros/primitives.mh
  // probably should be checked
  //
  //  enum SZINdir { BACKWARD= 0 ,  FORWARD = 1 }  ;

  
  if (Nd == 4)
  {
    c_1l = (Real)(5) / (Real)(8);
    c_3l = (Real)(-1) / (u0*u0*(Real)(16));
    c_5l = - c_3l / (u0*u0*(Real)(4));
    c_7l = - c_5l / (u0*u0*(Real)(6));
    c_Lepage = c_3l / (u0*u0);

                    
    for(mu=0; mu < Nd; ++mu)
      {
	uf[mu] = u[mu] * c_1l;
      
	for(nu=0; nu < Nd; ++nu) 
	  if(nu != mu)
	    {
	  tmp_0 = u[nu] * shift(u[mu], FORWARD, nu); 
	  tmp_2 = tmp_0 * shift(adj(u[nu]),  FORWARD, mu);
					
	  tmp_0 = u[nu] * shift(tmp_2,  FORWARD, nu);
	  tmp_1 = tmp_0 * shift(adj(u[nu]),  FORWARD, mu);

	  uf[mu] += tmp_1 * c_Lepage;

	  tmp_0 = u[mu] * shift(u[nu] , FORWARD, mu);
	  tmp_1 = shift(adj(u[nu]),  BACKWARD, nu) * shift(tmp_0, BACKWARD, nu);
			
	  tmp_2 += tmp_1;
	  uf[mu] += tmp_2 * c_3l;

	  tmp_0 = tmp_1 * shift(u[nu], FORWARD, mu);
	  tmp_1 = shift(adj(u[nu]),  BACKWARD, nu) * shift(tmp_0, BACKWARD, nu);
			
	  uf[mu] += tmp_1 * c_Lepage;
			
	  for(rho=0; rho < Nd; ++rho) if(rho != mu && rho != nu)
	  {
	    tmp_0 = u[rho] * shift(tmp_2, FORWARD, rho);
	    tmp_3 = tmp_0 * shift(adj(u[rho]), FORWARD, mu);

	    tmp_0 = tmp_2 * shift(u[rho], FORWARD, mu);
	    tmp_3 += shift(adj(u[rho]), BACKWARD, rho) * shift(tmp_0, BACKWARD, rho);
			   
	    uf[mu] += tmp_3 * c_5l;

	    for(sigma=0; sigma < Nd; ++sigma)
	    if(sigma != mu && sigma != nu && sigma != rho)
	    {
	      tmp_0 = u[sigma] * shift(tmp_3, FORWARD, sigma);
	      tmp_1 = tmp_0 * shift(adj(u[sigma]), FORWARD, mu);

	      tmp_0 = tmp_3 * shift(u[sigma], FORWARD, mu);
	      tmp_1 += shift(adj(u[sigma]), BACKWARD, sigma) * shift(tmp_0,  BACKWARD, sigma);

	      uf[mu] += tmp_1 * c_7l;
			     
	    }	  
	    }	     
	    }	   
            }	     
            }			    
                  
  else
    QDP_error_exit("Fat7_links not implemented for this dim", Nd);
  
  END_CODE();
}

}; // End Namespace Chroma

using namespace Chroma;

