// $Id: kyuqprop_io.h,v 1.3 2004-05-03 20:28:27 edwards Exp $
/*! \file
 *  \brief Read/write a KYU propagator
 */

#ifndef __kyuqprop_io_h__
#define __kyuqprop_io_h__

//! Initialize the Dirac to Degrand-Rossi spin transformation matrix
/*!
 * \ingroup io
 *
 * Initialize the similarity transformation matrix from 
 * Euclidean Dirac to Euclidean Degrand-Rossi basis
 *
 * \param U          spin matrix ( Modify )
 *
 * \returns The U in   Gamma_{Degrand-Rossi} = U Gamma_Dirac U^dag
 */    
void initDiracToDRMat(SpinMatrix& U);


//! Read a KYU propagator file
/*!
 * \ingroup io
 *
 * \param q          propagator ( Read )
 * \param file       path ( Read )
 */    
void readKYUQprop(LatticePropagator& q, const string& file);


#endif
