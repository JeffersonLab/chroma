// -*- C++ -*-
// $Id: hotst.h,v 1.2 2003-11-25 04:18:56 edwards Exp $
/*! \file
 *  \brief Set a gauge field from a sample of (almost) Haar measure
 */

#ifndef HOTST_INCLUDE
#define HOTST_INCLUDE

//! Set a gauge field from a sample of (almost) Haar measure
/*!
 * \ingroup gauge
 *
 * u = exp(A), where A = random traceless antihermitian matrix.
 *
 * Arguments:
 *
 *  \param u          Gauge field               (Modify)
 *
 *
 * NOTE!!!!: we are using a hack. The twelth_order exponentiation followed
 *   by a reunitarization should be replace by an exact exponentiation.
 *   So, the gauge field distriubution is NOT from Haar measure, but
 *   if someone can think of a really good reason why it should be Haar
 *   measure, then something can be done about it. 
 */
void HotSt(multi1d<LatticeColorMatrix>& u);

#endif
