// -*- C++ -*-
// $Id: phfctr.h,v 1.1 2003-12-16 22:18:48 edwards Exp $
/*! \file
 *  \brief  This routine is specific to Wilson fermions!
 */

#ifndef __setph_h__
#define __setph_h__

//!  Initialize phase factors
/*!
 * \ingroup gauge
 *
 * Setup the "phases" used by routine PHFCTR that multiplies the gauge
 * fields so as to handle boundary conditions.
 *
 * Arguments:
 *
 *  \param boundary       boundary conditions of the lattice       (Read)
 */

void setph(const multi1d<int>& boundary);


//! Multiply the gauge fields by the phase factors (-1)^X
/*!
 * \ingroup gauge
 *
 * Arguments:
 *
 *  \param u          Gauge field               (Modify)
 */
void phfctr(multi1d<LatticeColorMatrix>& u);


#endif
