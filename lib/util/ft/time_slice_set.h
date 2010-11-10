// -*- C++ -*-
/*! \file
 *  \brief Convenience for building time-slice subsets
 */

#ifndef __time_slice_set_h__
#define __time_slice_set_h__

#include "chromabase.h"

namespace Chroma 
{

  //! Builds time slice subsets
  /*!
   * \ingroup ft
   */
  class TimeSliceSet
  {
  public:
    //! Constructor about origin
    TimeSliceSet(int decay_dir);
    
    //! The set to be used in sumMulti
    const Set& getSet() const {return sft_set;}

    //! Number of subsets - length in decay direction
    int numSubsets() const {return sft_set.numSubsets();}

    //! Number of sites in each subset
    int numSites() const;

    //! Decay direction
    int getDir() const {return decay_dir;}

  private:
    TimeSliceSet() {} // hide default constructor

    int  decay_dir;
    Set  sft_set;
  };

}  // end namespace Chroma

#endif
