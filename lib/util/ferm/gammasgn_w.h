// -*- C++ -*-
// $Id: gammasgn_w.h,v 1.1 2004-06-11 20:31:25 edwards Exp $
/*! \file
 *  \brief Compute gamma matrix multiplication table factors
 */

#ifndef __gammasgn_h__
#define __gammasgn_h__

//! Return gamma matrix multiplication table factors
/*!
 * \ingroup ferm
 *
 * Initialize signs needed for  Gamma(n)*Gamma(m)=sgn(n,m)*Gamma(n ^ m)
 */

int gammaSgn(int n, int m);

#endif
