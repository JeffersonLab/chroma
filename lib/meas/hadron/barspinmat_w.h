// -*- C++ -*-
// $Id: barspinmat_w.h,v 3.4 2008-01-11 18:14:18 kostas Exp $
/*! \file
 *  \brief Baryon spin and projector matrices
 */

#ifndef __barspinmat_w_h__
#define __barspinmat_w_h__

#include "chromabase.h"

namespace Chroma 
{

  //! All the spin and projectors for 2 and 3-pts.
  /*! \ingroup hadron */
  namespace BaryonSpinMats
  { 
    //! NR = (1/2)* ( 1 + g_4 )
    SpinMatrix NR();

    //! NRnegPar = (1/2)* ( 1 - g_4 )
    SpinMatrix NRnegPar();

    //! C = Gamma(10)
    SpinMatrix C();

    //! C NR = (1/2)*C * ( 1 + g_4 )
    SpinMatrix CNR();

    //! C gamma_5 gamma_4 = - Gamma(13)
    SpinMatrix Cg5g4();

    //! C g_k = C gamma_k
    SpinMatrix Cgk(int k);

    //! C g_\mu = C gamma_\mu
    SpinMatrix Cgmu(int k);
    //! C g_\mu = C gamma_\mu
    SpinMatrix CgmuTrans(int k);

    //! C g4 g_k = C gamma_4 gamma_k
    SpinMatrix Cg4gk(int k);

    //! C g_k NR = C gamma_k (1/2)(1 + gamma_4)
    SpinMatrix CgkNR(int k);

    //! C g_5 = C gamma_5 = Gamma(5)
    SpinMatrix Cg5();

    //! C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
    SpinMatrix Cg5NR();

    //! C g_5 NR = (1/2)*C gamma_5 * ( 1 - g_4 )
    SpinMatrix Cg5NRnegPar();

    //! C gamma_- = Cgm = (C gamma_-)^T
    SpinMatrix Cgm();

    //! C gamma_4 gamma_- = Cg4m
    SpinMatrix Cg4m();

    //! C gamma_- NR = CgmNR = C gamma_- (1/2)(1 + gamma_4)
    SpinMatrix CgmNR();

   
    //! C gamma_+ = Cgp = (C gamma_+)^T
    SpinMatrix Cgp();

    //! C gamma_4 gamma_+ = Cg4p
    SpinMatrix Cg4p();

    //! C gamma_+ NR = CgpNR = C gamma_+ (1/2)(1 + gamma_4)
    SpinMatrix CgpNR();


    //! T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2
    SpinMatrix Tunpol();

    //! T = (1 + gamma_4) / 2 = (1 - Gamma(8)) / 2
    SpinMatrix TunpolNegPar();

    //! T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2
    SpinMatrix Tpol();

    //! T = \Sigma_3 (1 - gamma_4) / 2 = -i (-Gamma(3) + Gamma(11)) / 2
    SpinMatrix TpolNegPar();

    //! T = (1 + \Sigma_3)*(1 + gamma_4) / 2   = (1 + Gamma(8) - i G(3) - i G(11)) / 2
    SpinMatrix Tmixed();


    //! T = (1 - \Sigma_3)*(1 - gamma_4) / 2   = (1 - Gamma(8) - i G(3) + i G(11)) / 2
    // Need to flip the spin for time reversal
    SpinMatrix TmixedNegPar();
    
    //! T = (1 + i \gamma_5 \gamma_3 )/2 = 1/2 *( 1  - iG(11) )
    SpinMatrix TspinUp() ;
    //! T = (1 - i \gamma_5 \gamma_3 )/2 = 1/2 *( 1  + iG(11) )
    SpinMatrix TspinDown() ;
    //! T = i \gamma_5 (\gamma_1 + i \gamma_2 )
    SpinMatrix T_ig5XpiY();
    //! T = i \gamma_5 (\gamma_1 - i \gamma_2 )
    SpinMatrix T_ig5XmiY() ;

  }

}  // end namespace Chroma


#endif
