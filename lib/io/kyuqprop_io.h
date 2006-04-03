// $Id: kyuqprop_io.h,v 3.0 2006-04-03 04:58:55 edwards Exp $
/*! \file
 *  \brief Read/write a KYU propagator
 */

#ifndef __kyuqprop_io_h__
#define __kyuqprop_io_h__

namespace Chroma {

//! Read a KYU propagator file
/*!
 * \ingroup io
 *
 * \param q          propagator ( Read )
 * \param file       path ( Read )
 */    
void readKYUQprop(LatticePropagator& q, const string& file);

}  // end namespace Chroma

#endif
