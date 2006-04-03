// $Id: readszinqprop_w.h,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 *  \brief Read a SZIN propagator
 */

#ifndef __readszinqprop_h__
#define __readszinqprop_h__

namespace Chroma {

//! Read a SZIN propagator file. This is a simple memory dump readr.
/*!
 * \param xml        xml readr holding prop info ( Read )
 * \param q          propagator ( Read )
 * \param file       path ( Read )
 */    

void readSzinQprop(XMLReader& xml, LatticePropagator& q, const string& file);

//! Read a SZIN propagator file. This is a simple memory dump readr.
/*!
 * \ingroup io
 *
 * \param q          propagator ( Read )
 * \param file       path ( Read )
 * \param kappa      kappa value (Read)
 */    

void readSzinQprop(LatticePropagator& q, const string& file,
		   const Real& kappa);

}  // end namespace Chroma

#endif
