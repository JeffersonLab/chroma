// -*- C++ -*-
// $Id: taproj.h,v 2.0 2005-09-25 21:04:45 edwards Exp $
/*! \file
 *  \brief Take the traceless antihermitian projection of a color matrix
 */

#ifndef __taproj_h__
#define __taproj_h__

namespace Chroma {

//! Take the traceless antihermitian projection of a color matrix
/*!
 * \ingroup gauge
 *
 *  a = (1/2)[a - a_dag] - Tr[(1/2)*(a - a_dag)]/Nc
 *
 * that is the anti-hermitian traceless part of a 
 *
 * Arguments:
 *
 *  \param a        LatticeColorMatrix          (Modify)
 */

void taproj(LatticeColorMatrix& a);

}  // end namespace Chroma

#endif
