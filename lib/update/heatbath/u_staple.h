// -*- C++ -*-
// $Id: u_staple.h,v 3.0 2006-04-03 04:59:07 edwards Exp $
/*! \file
 *  \brief calculate the sum of all staples for a particular u_mu link
 */

#ifndef __u_staple__
#define __u_staple__

#include "update/heatbath/hb_params.h"

namespace Chroma 
{

  //! Compute staple from Wilson gauge action
  /*!
   * \ingroup heatbath
   *
   * \param u                    link field
   * \param mu                   direction of the link
   * \param u_mu_staple          staple attached to the link u
   * \param sub                  subset for updating
   * \param hpb                  container of HB parameters
   */
  void u_staple(LatticeColorMatrix& u_mu_staple,
                const multi1d<LatticeColorMatrix>& u,
                int mu,
                const OrderedSubset& sub,
                const HBParams& hbp);

}  // end namespace Chroma

#endif
