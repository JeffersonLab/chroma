/*! \file
 *  \brief Convenience for building time-slice subsets
 */

#include "util/ft/time_slice_set.h"
#include "util/ft/single_phase.h"
#include "qdp_util.h"                 // part of QDP++, for crtesn()

namespace Chroma 
{

  // Anonymous namespace
  namespace
  {
    //! Function object used for constructing the time-slice set
    class TimeSliceFunc : public SetFunc
    {
    public:
      TimeSliceFunc(int dir) : dir_decay(dir) {}

      int operator() (const multi1d<int>& coordinate) const
      {
	if ((dir_decay<0)||(dir_decay>=Nd)) {
	  return 0 ;
	} else {
	  return coordinate[dir_decay] ;
	}
      }

      int numSubsets() const
      {
	if ((dir_decay<0)||(dir_decay>=Nd)) {
	  return 1 ;
	} else {
	  return Layout::lattSize()[dir_decay] ;
	}
      }

    private:
      TimeSliceFunc() {}  // hide default constructor

      int dir_decay;
    };

  } // end anonymous namespace



  int
  TimeSliceSet::numSites() const
  {
    int vol = 1;

    if ((decay_dir<0)||(decay_dir>=Nd))
      vol = Layout::vol();
    else 
    {
      for(int m=0; m < Nd; ++m)
	vol *= Layout::lattSize()[m];
    }

    return vol;
  }


  TimeSliceSet::TimeSliceSet(int j_decay) : decay_dir(j_decay)
  {
    sft_set.make(TimeSliceFunc(j_decay));
  }


}  // end namespace Chroma
