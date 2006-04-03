// -*- C++ -*-
// $Id: paulitodr.h,v 3.0 2006-04-03 04:59:11 edwards Exp $
/*! \file
 *  \brief Basis rotation matrix from Pauli-Schwinger (Euclidean Sakurai) to Degrand-Rossi (and reverse)
 */

#ifndef __paulitodr_h__
#define __paulitodr_h__

namespace Chroma {

//! The Pauli-Schwinger (Euclidean Sakurai) to Degrand-Rossi spin transformation matrix
/*!
 * \ingroup ferm
 *
 * Return the similarity transformation matrix from 
 * Euclidean Pauli-Schwinger to Euclidean Degrand-Rossi basis
 *
 * \returns The U in   Gamma_{Degrand-Rossi} = U Gamma_PS U^dag
 */

SpinMatrixD PauliToDRMat();

}  // end namespace Chroma

#endif
