// $Id: writeszinqprop_w.h,v 1.1 2003-09-25 22:20:28 edwards Exp $

#ifndef __writeszinqprop_h__
#define __writeszinqprop_h__

/*! \file
 *  \brief Write out a SZIN propagator
 */

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


#endif
