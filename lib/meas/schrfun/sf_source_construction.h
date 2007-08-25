// -*- C++ -*-
// $Id: sf_source_construction.h,v 1.1 2007-08-25 04:07:41 edwards Exp $

/*! @file
 * @brief Source construction within Schroedinger Functional
 */

#ifndef __sf_source_construction_h__
#define __sf_source_construction_h__

#include "chromabase.h"

namespace Chroma
{
  //! Parameters used for SF source construction
  /*! @ingroup sources
   *
   * Supports creation of quark sources
   */
  struct SFSourceConstParams_t
  {
    int t0;             /*!< Time slice source location */
    int decay_dir;      /*!< Decay direction */
    int color_source;   /*!< Color source */
    int spin_source;    /*!< Color source */
  };


  //! Base class for quark source construction
  /*! @ingroup sources
   *
   * Supports creation of quark sources
   */
  template<typename T>
  class SFSourceConstruction
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SFSourceConstruction() {}

    //! Construct the source on a fixed time-slice and decay direction
    virtual T operator()(const multi1d<LatticeColorMatrix>& u, 
			 int t0, int decay_dir) const = 0;
  };

}


#endif
