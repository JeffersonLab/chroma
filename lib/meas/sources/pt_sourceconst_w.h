// -*- C++ -*-
// $Id: pt_sourceconst_w.h,v 1.1 2005-10-28 21:06:41 edwards Exp $
/*! \file
 *  \brief Point source construction
 */

#ifndef __pt_sourceconst_w_h__
#define __pt_sourceconst_w_h__

#include "meas/sources/source_construction.h"
#include "meas/sources/ptsrc_params.h"

namespace Chroma
{

  //! Name and registration
  namespace PropPointSourceConstEnv
  {
    extern const std::string name;
    extern const bool registered;
    //! Name to be used
  }
  

  //! Point source construction
  /*! @ingroup sources
   *
   * Create a point propagator source
   */
  class PropPointSourceConst : public SourceConstruction<LatticePropagator>
  {
  public:
    //! Full constructor
    PropPointSourceConst(const PointSourceConstParams& p) : params(p) {}

    //! Construct the source
    LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u) const;

  private:
    //! Hide partial constructor
    PropPointSourceConst() {}

  private:
    PointSourceConstParams  params;   /*!< source params */
  };

}  // end namespace Chroma


#endif
