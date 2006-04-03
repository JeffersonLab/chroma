// $Id: seqpiontest_w.cc,v 3.0 2006-04-03 04:59:00 edwards Exp $
/*! \file
 *  \brief Test a pion sequential source
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/seqpiontest_w.h"

namespace Chroma {

//! Test a pion sequential source.
/*!
 * \ingroup hadron
 *
 *  For the case of a pion, we have evaluated as the sequential source
 *
 *  H(y, 0; tx, p) = \sum exp{ip.x} U(y,x) \gamma_5 D(x,0) \gamma_5
 *
 *  Thus we can see that 
 *
 *  Tr H(0,0; tx, p) = \sum_x exp{ip.x} Tr[ U^dagger(x,0) D(x,0)]
 * 
 *  which is just the conjugate of the pion correlator at momentum
 *  p and timslice tx
 *
 *  ARGUMENTS:
 *
 *  \param pion_src          the sequential propagator evaluated at t_source (Write)
 *  \param seq_quark_prop    the sequential propagator ( Read )
 *  \param t_source          the coordinates of the source ( Read )
 */

void seqPionTest(Complex& pion_src,
		 const LatticePropagator& seq_quark_prop,
		 const multi1d<int>& t_source)
{
  START_CODE();

  /*
   *  We now pick out the particular component of the quark propagator
   *  at the source, and take the trace
   */

  pion_src = trace(peekSite(seq_quark_prop, t_source));

  END_CODE();
}
  
}  // end namespace Chroma
