/*! \file
 * \brief LatticeColorVector time-slice IO cache
 *
 * LatticeColorVector time-slice IO cache
 */

#include "util/ferm/timeslice_io_cache.h"

namespace Chroma 
{  
  //----------------------------------------------------------------------------
  // Constructor
  TimeSliceIOCache::TimeSliceIOCache(QDP::MapObjectDisk< KeyTimeSliceColorVec_t,TimeSliceIO<LatticeColorVector> >& eigen_source_)
    : eigen_source(eigen_source_)
  {
    const int Lt = Layout::lattSize()[Nd-1];

    // Figure out how many vectors are in the source
    // We know time slice 0 has to be a part of the sources
    num_vecs = 0;
    while(1)
    {
      KeyTimeSliceColorVec_t key;
      key.t_slice  = 0;
      key.colorvec = num_vecs;

      if (! eigen_source.exist(key)) {break;}

      ++num_vecs;
    }

    if (num_vecs == 0)
    {
      QDPIO::cerr << __func__ << ": this is bad - did not find any eigenvectors in eigen_source\n";
      QDP_abort(1);
    }
    else
    {
      QDPIO::cout << __func__ << ": found in eigenvector source num_vecs= " << num_vecs << std::endl;
    }

    eigen_cache.resize(num_vecs);
    cache_marker.resize(Lt,num_vecs);

    for(int n=0; n < num_vecs; ++n)
    {
      eigen_cache[n] = zero;
      
      for(int t=0; t < Lt; ++t)
	cache_marker(t,n) = false;
    }
  }


  // Get a vector
  LatticeColorVector& TimeSliceIOCache::getVec(int colorvec)
  {
    return eigen_cache[colorvec];
  }

  // Get a vector
  LatticeColorVector& TimeSliceIOCache::getVec(int t_actual, int colorvec)
  {
    // If not in cache, then retrieve
    if (! cache_marker(t_actual,colorvec))
    {
      KeyTimeSliceColorVec_t key_vec;
      key_vec.t_slice  = t_actual;
      key_vec.colorvec = colorvec;

      TimeSliceIO<LatticeColorVector> time_slice_io(eigen_cache[colorvec], t_actual);

      eigen_source.get(key_vec, time_slice_io);
      cache_marker(t_actual,colorvec) = true;
    }

    return eigen_cache[colorvec];
  }

} // namespace Chroma
