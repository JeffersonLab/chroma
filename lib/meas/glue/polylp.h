// -*- C++ -*-
// $Id: polylp.h,v 1.2 2003-10-10 03:46:46 edwards Exp $
/*! \file
 *  \brief Calculate the global normalized sum of the Polyakov loop
 */

#ifndef __polylp_h__
#define __polylp_h__

//! Return the value of the average plaquette normalized to 1
/*!
 * \ingroup glue
 *
 * \param u          gauge field (Read)
 * \param poly_loop  Polyakov loop average in direction mu (Write) 
 * \param mu         direction of Polyakov loop (Read)
 */

void polylp(const multi1d<LatticeColorMatrix>& u, DComplex& poly_loop, int mu);

#endif
