// $Id: readszinqprop_w.h,v 1.4 2003-06-20 20:48:12 dgr Exp $

#ifndef __readszinqprop_h__
#define __readszinqprop_h__

/*! \file
 *  \brief Read in a configuration written by SZIN up to configuration version 7.
 */

//! Read a SZIN propagator file. This is a simple memory dump reader.
/*!
 * \param q          propagator ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readSzinQprop(LatticePropagator& q, const string& file);

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
