// $Id: sink_smear2_w.cc,v 1.1 2003-03-07 05:42:01 edwards Exp $
/*! \file
 *  \brief Control routine for types of propagator smearing
 */

#include "chromabase.h"
#include "primitives.h"
#include "meas/smear/gaus_smear.h"

using namespace QDP;

//! "Smear" the quark propagator at the sink by a covariant Gaussian
/*!
 * This routine is specific to Wilson fermions!
 *
 * Arguments:
 *
 *  \param u                   gauge field ( Read )
 *  \param quark_propagator    quark propagator ( Modify )
 *  \param wvf_type            wave function type: Gaussian or exponential ( Read )
 *  \param wvf_param           wvf_param of "shell" wave function ( Read )
 *  \param WvfIntPar           number of iterations to approximate Gaussian
 *                             or terminate CG inversion for Wuppertal smearing ( Read )
 *  \param j_decay             direction of decay ( Read ) 
 */

void sink_smear2(const multi1d<LatticeColorMatrix>& u, 
		 LatticePropagator& quark_propagator, 
		 int wvf_type, const Real& wvf_param, int WvfIntPar, int j_decay)
{
  Real RsdCG;

  START_CODE("sink_smear2");

  switch (wvf_type)
  {
  case OPTION_GAUGE_INV_GAUSSIAN_WVF:
    gausSmear(u, quark_propagator, wvf_param, WvfIntPar, j_decay);
    break;

#if 0
  case OPTION_WUPPERTAL_WVF:
    RsdCG = fuzz;
    wupp_smear(u, quark_propagator, wvf_param, WvfIntPar, j_decay, RsdCG);
    break;
#endif

  default:
    QDP_error_exit("Unknown gauge invariant wave function", wvf_type);
  }
    
  END_CODE("sink_smear2");
}
