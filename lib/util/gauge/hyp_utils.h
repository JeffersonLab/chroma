// -*- C++ -*-
/*! \file
 *  \brief Hyp utilities
 */

#ifndef HYP_UTILS_H
#define HYP_UTILS_H

#include "chromabase.h"


namespace Chroma 
{

  /*!
   * Hyping
   *
   * \ingroup gauge
   *
   * @{
   */

  /*! \ingroup gauge */
  namespace HypLinkTimings { 
    double getForceTime();
    double getSmearingTime();
    double getFunctionsTime();
  }

  /*! \ingroup gauge */
  namespace Hyping 
  {
    //! Do the smearing from level i to level i+1
    void smear_links(const multi1d<LatticeColorMatrix>& current, 
                     multi1d<LatticeColorMatrix>& next,
                     multi1d<LatticeColorMatrix>& Omega1,
                     multi1d<LatticeColorMatrix>& Omega2,
                     multi1d<LatticeColorMatrix>& Omega3,
                     multi1d<LatticeColorMatrix>& QPowHalf1,
                     multi1d<LatticeColorMatrix>& QPowHalf2,
                     multi1d<LatticeColorMatrix>& QPowHalf3,
                     const multi1d<bool>& smear_in_this_dirP,
                     const Real alpha1,
                     const Real alpha2,
                     const Real alpha3,
                     const int BlkMax,
                     const Real BlkAccu);

    /*! \ingroup gauge */
    void hyp_lv1_links(const multi1d<LatticeColorMatrix>& u, 
                       multi1d<LatticeColorMatrix>& u_lv1,
                       multi1d<LatticeColorMatrix>& Omega,
                       multi1d<LatticeColorMatrix>& QPowHalf,
                       const multi1d<bool>& smear_in_this_dirP,
                       const Real alpha1,
                       const Real alpha2,
                       const Real alpha3,
                       const int BlkMax,
                       const Real BlkAccu);

    /*! \ingroup gauge */
    void hyp_lv2_links(const multi1d<LatticeColorMatrix>& u, 
                       multi1d<LatticeColorMatrix>& u_lv1,
                       multi1d<LatticeColorMatrix>& u_lv2,
                       multi1d<LatticeColorMatrix>& Omega,
                       multi1d<LatticeColorMatrix>& QPowHalf,
                       const multi1d<bool>& smear_in_this_dirP,
                       const Real alpha1,
                       const Real alpha2,
                       const Real alpha3,
                       const int BlkMax,
                       const Real BlkAccu);

    /*! \ingroup gauge */
    void hyp_lv3_links(const multi1d<LatticeColorMatrix>& u, 
                       multi1d<LatticeColorMatrix>& u_lv3,
                       multi1d<LatticeColorMatrix>& u_hyp,
                       multi1d<LatticeColorMatrix>& Omega,
                       multi1d<LatticeColorMatrix>& QPowHalf,
                       const multi1d<bool>& smear_in_this_dirP,
                       const Real alpha1,
                       const Real alpha2,
                       const Real alpha3,
                       const int BlkMax,
                       const Real BlkAccu);

    //! Do the force recursion from level i+1, to level i
    void deriv_recurse(multi1d<LatticeColorMatrix>&  F,
		       const multi1d<bool>& smear_in_this_dirP,
                       const Real alpha1,
                       const Real alpha2,
                       const Real alpha3,
                       const int hyp_qr_max_iter,
                       const Real hyp_qr_tol,
                       const multi1d<LatticeColorMatrix>& u);
    
    //! Compute Upper Hessenberg reduction
    void upper_hessenberg(const LatticeColorMatrix &U, 
                          LatticeColorMatrix &UH);
    
    //! Compute QR reduction from Upper Hessenberg reduction
    void qr_from_upper_hess(LatticeColorMatrix &UH,
                            const Real hyp_qr_tol,
                            const int hyp_qr_max_iter);

    //! Compute vandermonde inversion
    void solve_vandermonde(LatticeColorMatrix &UT,
                           multi1d<LatticeComplex>& f);
    
  }
  
  /*! @} */   // end of group gauge
}

#endif
