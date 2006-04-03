// $Id: writeszinqprop_w.h,v 3.0 2006-04-03 04:58:56 edwards Exp $

#ifndef __writeszinqprop_h__
#define __writeszinqprop_h__

/*! \file
 *  \brief Write out a SZIN propagator
 */

namespace Chroma {


//! Write a SZIN propagator file. This is a simple memory dump writer.
/*!
 * \ingroup io
 *
 * \param q          propagator ( Read )
 * \param file       path ( Read )
 * \param kappa      kappa value (Read)
 */    

void writeSzinQprop(const LatticePropagator& q, const string& file,
		    const Real& kappa);

}  // end namespace Chroma

#endif
