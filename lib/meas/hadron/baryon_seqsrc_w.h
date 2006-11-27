// -*- C++ -*-
// $Id: baryon_seqsrc_w.h,v 3.1 2006-11-27 04:33:35 edwards Exp $
/*! \file
 *  \brief Construct baryon sequential sources.
 */

#ifndef __baryon_seqsrc_w_h__
#define __baryon_seqsrc_w_h__

#include "meas/hadron/hadron_seqsource.h"

namespace Chroma 
{

  //! Baryon-Baryon seqsources have a time order phase
  /*! @ingroup hadron */
  class BaryonSeqSourceBase : public HadronSeqSource<LatticePropagator>
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~BaryonSeqSourceBase() {}

    //! Construct the source
    virtual LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
					 const multi1d<ForwardProp_t>& forward_headers,
					 const multi1d<LatticePropagator>& forward_props) = 0;
    
  protected:
    //! Combine projection with time-ordering
    virtual LatticePropagator projectBaryon(const LatticePropagator& src_prop_tmp,
					    const multi1d<ForwardProp_t>& forward_headers);

    //! Time-ordering phase of source and sink hadron states
    virtual Complex timeOrder() const;

    //! Convenience function to yank the boundary condition from the forward prop headers
    virtual void setBC(const multi1d<ForwardProp_t>& forward_headers);

    //! Set bc
    virtual multi1d<int>& getBC() = 0;

    //! Get bc
    virtual const multi1d<int>& getBC() const = 0;
  };

}  // end namespace Chroma


#endif
