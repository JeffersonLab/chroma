// -*- C++ -*-
// $Id: barseqsrc_w.h,v 2.0 2005-09-25 21:04:35 edwards Exp $
/*! \file
 *  \brief Construct baryon sequential sources.
 */

#ifndef __barseqsrc_w_h__
#define __barseqsrc_w_h__

namespace Chroma 
{

  LatticePropagator barNuclUUnpol(const LatticePropagator& quark_propagator_1, 
				  const LatticePropagator& quark_propagator_2,
				  const LatticePropagator& quark_propagator_3);

  LatticePropagator barNuclDUnpol(const LatticePropagator& quark_propagator_1, 
				  const LatticePropagator& quark_propagator_2,
				  const LatticePropagator& quark_propagator_3);

  LatticePropagator barNuclUPol(const LatticePropagator& quark_propagator_1, 
				const LatticePropagator& quark_propagator_2,
				const LatticePropagator& quark_propagator_3);

  LatticePropagator barNuclDPol(const LatticePropagator& quark_propagator_1, 
				const LatticePropagator& quark_propagator_2,
				const LatticePropagator& quark_propagator_3);

  LatticePropagator barDeltaUUnpol(const LatticePropagator& quark_propagator_1, 
				   const LatticePropagator& quark_propagator_2,
				   const LatticePropagator& quark_propagator_3);

  LatticePropagator barDeltaDUnpol(const LatticePropagator& quark_propagator_1, 
				   const LatticePropagator& quark_propagator_2,
				   const LatticePropagator& quark_propagator_3);

  LatticePropagator barNuclUUnpolNR(const LatticePropagator& quark_propagator_1, 
				    const LatticePropagator& quark_propagator_2,
				    const LatticePropagator& quark_propagator_3);

  LatticePropagator barNuclDUnpolNR(const LatticePropagator& quark_propagator_1, 
				    const LatticePropagator& quark_propagator_2,
				    const LatticePropagator& quark_propagator_3);

  LatticePropagator barNuclUPolNR(const LatticePropagator& quark_propagator_1, 
				  const LatticePropagator& quark_propagator_2,
				  const LatticePropagator& quark_propagator_3);

  LatticePropagator barNuclDPolNR(const LatticePropagator& quark_propagator_1, 
				  const LatticePropagator& quark_propagator_2,
				  const LatticePropagator& quark_propagator_3);

  LatticePropagator barNuclUMixedNR(const LatticePropagator& quark_propagator_1, 
				    const LatticePropagator& quark_propagator_2,
				    const LatticePropagator& quark_propagator_3);

  LatticePropagator barNuclDMixedNR(const LatticePropagator& quark_propagator_1, 
				    const LatticePropagator& quark_propagator_2,
				    const LatticePropagator& quark_propagator_3);

  //! Patch for the quarkContract12 piece in NuclUMixedNR and NuclDMixedNR
  LatticePropagator barNuclPatchMixedNR(const LatticePropagator& quark_propagator_1, 
					const LatticePropagator& quark_propagator_2,
					const LatticePropagator& quark_propagator_3);

}  // end namespace Chroma


#endif
