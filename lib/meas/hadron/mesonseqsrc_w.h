// -*- C++ -*-
// $Id: mesonseqsrc_w.h,v 1.4 2005-03-07 02:55:20 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#ifndef __mesonseqsrc_w_h__
#define __mesonseqsrc_w_h__

namespace Chroma 
{

  LatticePropagator mesPionSeqSrc(const LatticePropagator& quark_propagator_1, 
				  const LatticePropagator& quark_propagator_2,
				  const LatticePropagator& quark_propagator_3);

}  // end namespace Chroma

#endif
