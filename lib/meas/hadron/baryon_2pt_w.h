// -*- C++ -*-
// $Id: baryon_2pt_w.h,v 1.1 2007-05-09 17:19:44 edwards Exp $
/*! \file
 *  \brief Construct baryon 2pt correlators.
 */

#ifndef __baryon_2pt_w_h__
#define __baryon_2pt_w_h__

#include "meas/hadron/hadron_2pt.h"

namespace Chroma 
{

  //! Baryon-Baryon seqsources have a time order phase
  /*! @ingroup hadron */
  class Baryon2PtBase : public HadronCorrelator
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~Baryon2PtBase() {}

    //! Construct the correlators
    virtual multi1d<Hadron2PtContraction_t> operator()(const multi1d<LatticeColorMatrix>& u) = 0;
    
  protected:
    //! Convenience function to yank the boundary condition from the forward prop headers
    virtual void setBC(const multi1d<ForwardProp_t>& forward_headers);

    //! Set bc
    virtual multi1d<int>& getBC() = 0;

    //! Get bc
    virtual const multi1d<int>& getBC() const = 0;
  };

}  // end namespace Chroma


#endif
