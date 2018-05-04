// -*- C++ -*-
/*! \file
 *  \brief Compute gamma matrix multiplication table factors
 */

#ifndef __gammasgn_w_h__
#define __gammasgn_w_h__

namespace Chroma 
{

  //! Return gamma matrix multiplication table factors
  /*!
   * \ingroup ferm
   *
   * Initialize signs needed for  Gamma(n)*Gamma(m)=sgn(n,m)*Gamma(n ^ m)
   */

  int gammaSgn(int n, int m);

}  // end namespace Chroma

#endif
