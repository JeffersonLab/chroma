// -*- C++ -*-
// $Id: diractodr.h,v 1.1 2004-05-14 00:23:41 edwards Exp $
/*! \file
 *  \brief Basis rotation matrix from Dirac to Degrand-Rossi (and reverse)
 */

#ifndef __diractodr_h__
#define __diractodr_h__

//! The Dirac to Degrand-Rossi spin transformation matrix
/*!
 * \ingroup ferm
 *
 * Return the similarity transformation matrix from 
 * Euclidean Dirac to Euclidean Degrand-Rossi basis
 *
 * \returns The U in   Gamma_{Degrand-Rossi} = U Gamma_Dirac U^dag
 */

SpinMatrix DiracToDRMat();

#endif
