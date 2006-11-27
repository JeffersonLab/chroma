// -*- C++ -*-
// $Id: meson_seqsrc_w.h,v 3.1 2006-11-27 04:33:35 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#ifndef __meson_seqsrc_w_h__
#define __meson_seqsrc_w_h__

#include "meas/hadron/hadron_seqsource.h"

namespace Chroma 
{

  //! Base class for meson sequential source construction
  /*! @ingroup hadron
   */
  class MesonSeqSourceBase : public HadronSeqSource<LatticePropagator>
  {
  public:
    //! Default destructor
    virtual ~MesonSeqSourceBase() {}
      
    //! Construct the source
    virtual LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
					 const multi1d<ForwardProp_t>& forward_headers,
					 const multi1d<LatticePropagator>& forward_props) = 0;
    
  };

}  // end namespace Chroma

#endif
