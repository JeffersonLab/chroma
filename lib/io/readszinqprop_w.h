// $Id: readszinqprop_w.h,v 1.8 2004-04-14 20:03:13 edwards Exp $
/*! \file
 *  \brief Read a SZIN propagator
 */

#ifndef __readszinqprop_h__
#define __readszinqprop_h__

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


#endif
