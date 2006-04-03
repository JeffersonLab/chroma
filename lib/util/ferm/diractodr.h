// -*- C++ -*-
// $Id: diractodr.h,v 3.0 2006-04-03 04:59:11 edwards Exp $
/*! \file
 *  \brief Basis rotation matrix from Dirac to Degrand-Rossi (and reverse)
 */

#ifndef __diractodr_h__
#define __diractodr_h__

namespace Chroma {

//! The Dirac to Degrand-Rossi spin transformation matrix
/*!
 * \ingroup ferm
 *
 * Return the similarity transformation matrix from 
 * Euclidean Dirac to Euclidean Degrand-Rossi basis
 *
 * \returns The U in   Gamma_{Degrand-Rossi} = U Gamma_Dirac U^dag
 */

SpinMatrixD DiracToDRMat();

}  // end namespace Chroma

#endif
