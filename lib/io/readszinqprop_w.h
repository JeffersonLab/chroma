// $Id: readszinqprop_w.h,v 1.3 2003-04-30 21:19:33 edwards Exp $

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

#endif
