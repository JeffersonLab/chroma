// -*- C++ -*-
// $Id: gammasgn_w.h,v 2.0 2005-09-25 21:04:43 edwards Exp $
/*! \file
 *  \brief Compute gamma matrix multiplication table factors
 */

#ifndef __gammasgn_h__
#define __gammasgn_h__

namespace Chroma {

//! Return gamma matrix multiplication table factors
/*!
 * \ingroup ferm
 *
 * Initialize signs needed for  Gamma(n)*Gamma(m)=sgn(n,m)*Gamma(n ^ m)
 */

int gammaSgn(int n, int m);

}  // end namespace Chroma

#endif
