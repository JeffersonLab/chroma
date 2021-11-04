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
                       const int BlkMax,
                       const Real BlkAccu,
                       const multi1d<LatticeColorMatrix>& u);
    
  }

  /*! @} */   // end of group gauge
}

#endif
