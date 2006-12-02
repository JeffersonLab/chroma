// -*- C++ -*-
// $Id: gammasgn_w.h,v 3.1 2006-12-02 18:13:37 edwards Exp $
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
