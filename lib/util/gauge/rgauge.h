// -*- C++ -*-
// $Id: rgauge.h,v 1.2 2003-10-08 18:54:37 edwards Exp $

#ifndef RGAUGE_INCLUDE
#define RGAUGE_INCLUDE

//! Do a random gauge transformation on the u fields
/*!
 * \ingroup gauge
 *
 * Arguments:
 *
 *  \param u          Gauge field                   (Modify)
 *  \param g          Gauge transformation matrices (Write)
 */

void rgauge(multi1d<LatticeColorMatrix>& u, LatticeColorMatrix& g);


//! Do a random gauge transformation on the u fields
/*!
 * \ingroup gauge
 *
 * Convenience function: does not return gauge transformation matrices
 *
 * Arguments:
 *
 *  \param u          Gauge field                   (Modify)
 */

void rgauge(multi1d<LatticeColorMatrix>& u);

#endif
