// -*- C++ -*-
// $Id: sink_smear2_w.h,v 1.7 2004-01-05 21:47:20 edwards Exp $
/*! \file
 *  \brief Control routine for types of propagator smearing
 */

#ifndef __sink_smear2_h__
#define __sink_smear2_h__

enum WvfType {
  WVF_TYPE_GAUGE_INV_GAUSSIAN,
  WVF_TYPE_WUPPERTAL,
  WVF_TYPE_UNKNOWN
};


//! "Smear" the quark propagator at the sink by a covariant Gaussian
/*!
 * \ingroup smear
 *
 * This routine is specific to Wilson fermions!
 *
 * Arguments:
 *
 *  \param u                   gauge field ( Read )
 *  \param quark_propagator    quark propagator ( Modify )
 *  \param wvf_type            wave function kind: Gaussian or exponential
 *                             ( Read )
 *  \param wvf_param           wvf_param of "shell" wave function ( Read )
 *  \param WvfIntPar           number of iterations to approximate Gaussian
 *                             or terminate CG inversion for Wuppertal smearing
 *                             ( Read )
 *  \param j_decay             direction of decay ( Read ) 
 */

void
sink_smear2(const multi1d<LatticeColorMatrix>& u,
            LatticePropagator& quark_propagator, WvfType wvf_type,
            const Real& wvf_param, int WvfIntPar, int j_decay);

#endif
