// $Id: readszinqprop_w.h,v 1.5 2003-08-27 22:08:41 edwards Exp $

#ifndef __readszinqprop_h__
#define __readszinqprop_h__

/*! \file
 *  \brief Read in a configuration written by SZIN up to configuration version 7.
 */

//! Read a SZIN propagator file. This is a simple memory dump reader.
/*!
 * \param xml        xml reader holding prop info ( Modify )
 * \param q          propagator ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readSzinQprop(XMLReader& xml, LatticePropagator& q, const string& file);

//! Write a SZIN propagator file. This is a simple memory dump writer.
/*!
 * \ingroup io
 *
 * \param q          propagator ( Read )
 * \param file       path ( Read )
 * \param kappa      kappa value (Read)
 */    

void writeSzinQprop(const LatticePropagator& q, const string& file,
		    const Real kappa);


#endif
