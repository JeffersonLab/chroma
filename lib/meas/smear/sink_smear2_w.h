// -*- C++ -*-
// $Id: sink_smear2_w.h,v 1.10 2004-01-13 03:58:32 edwards Exp $
/*! \file
 *  \brief Control routine for types of propagator smearing
 */

#ifndef __sink_smear2_h__
#define __sink_smear2_h__

#include "meas/smear/wvfkind.h"

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
 *  \param Wvf_kind            wave function kind: Gaussian or exponential
 *                             ( Read )
 *  \param wvf_param           wvf_param of "shell" wave function ( Read )
 *  \param WvfIntPar           number of iterations to approximate Gaussian
 *                             or terminate CG inversion for Wuppertal smearing
 *                             ( Read )
 *  \param j_decay             direction of decay ( Read ) 
 */

void
sink_smear2(const multi1d<LatticeColorMatrix>& u,
            LatticePropagator& quark_propagator, WvfKind Wvf_kind,
            const Real& wvf_param, int WvfIntPar, int j_decay);

#endif
