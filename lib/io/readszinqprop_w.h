// $Id: readszinqprop_w.h,v 1.7 2003-10-10 03:46:46 edwards Exp $

#ifndef __readszinqprop_h__
#define __readszinqprop_h__

/*! \file
 *  \brief Read out a SZIN propagator
 */

//! Read a SZIN propagator file. This is a simple memory dump readr.
/*!
 * \param xml        xml readr holding prop info ( Read )
 * \param q          propagator ( Read )
 * \param cfg_file   path ( Read )
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
