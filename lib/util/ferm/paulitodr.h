// -*- C++ -*-
// $Id: paulitodr.h,v 1.1 2004-05-25 16:17:01 edwards Exp $
/*! \file
 *  \brief Basis rotation matrix from Pauli-Schwinger (Euclidean Sakurai) to Degrand-Rossi (and reverse)
 */

#ifndef __paulitodr_h__
#define __paulitodr_h__

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

#endif
