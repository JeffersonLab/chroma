// $Id: kyuqprop_io.h,v 1.4 2004-05-14 00:21:46 edwards Exp $
/*! \file
 *  \brief Read/write a KYU propagator
 */

#ifndef __kyuqprop_io_h__
#define __kyuqprop_io_h__

//! Read a KYU propagator file
/*!
 * \ingroup io
 *
 * \param q          propagator ( Read )
 * \param file       path ( Read )
 */    
void readKYUQprop(LatticePropagator& q, const string& file);


#endif
