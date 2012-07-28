// -*- C++ -*-
/*! \file
 * \brief LatticeColorVector time-slice IO cache
 *
 * LatticeColorVector time-slice IO cache
 */

#ifndef __timeslice_io_cache_h__
#define __timeslice_io_cache_h__

#include "chromabase.h"
#include "qdp_map_obj_disk.h"
#include "util/ferm/key_timeslice_colorvec.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //! Cache for holding time slice eigenvectors
  class TimeSliceIOCache
  {
  public:
    //! Constructor
    TimeSliceIOCache(QDP::MapObjectDisk<KeyTimeSliceColorVec_t,LatticeColorVector>& eigen_source_);

    //! Virtual destructor
    virtual ~TimeSliceIOCache() {}

    //! Get number of vectors
    virtual int getNumVecs() const {return num_vecs;}

    //! Get the whole vector
    virtual LatticeColorVector& getVec(int colorvec);

    //! Get a vector
    virtual LatticeColorVector& getVec(int t_actual, int colorvec);

  private:
    // Arguments
    QDP::MapObjectDisk< KeyTimeSliceColorVec_t,TimeSliceIO<LatticeColorVector> >& eigen_source;

    // Local
    multi1d<LatticeColorVector>  eigen_cache;
    multi2d<bool>                cache_marker;
    int                          num_vecs;
  };

}

#endif
