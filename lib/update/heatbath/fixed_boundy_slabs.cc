/*! \file
 *  \brief Convenience for building time-slab subsets
 */

#include "fixed_boundy_slabs.h"
#include "qdp_util.h"                 // part of QDP++, for crtesn()

namespace Chroma 
{

  // Anonymous namespace
  namespace
  {
    //! Function object used for constructing the time-slab set
    class FixedBoundaryFunc : public SetFunc
    {
    public:
      FixedBoundaryFunc(int dir,multi1d<int> ti, multi1d<int> th,int N) : dir_decay(dir),t0(ti),thick(th),Nt(N) {}

      int operator() (const multi1d<int>& coordinate) const
      {

	if ((dir_decay<0)||(dir_decay>=Nd)) {
	  return 0 ;
	} else {
	  for(int k(0);k<t0.size();k++){
	      int t_eff = (coordinate[dir_decay] - t0[k] +  Nt) % Nt;
	      if(t_eff<abs(thick[k]))
		return 1 ;
	  }
	}
	return 0 ;
      }

      int numSubsets() const
      {
	if ((dir_decay<0)||(dir_decay>=Nd)) {
	  return 1 ;
	} else {
	  return 2 ;
	}
      }

    private:
      FixedBoundaryFunc() {}  // hide default constructor

      int dir_decay;
      multi1d<int> t0;
      multi1d<int> thick;
      int Nt ;
    };


  } // end anonymous namespace

  FixedBoundarySlabSet::FixedBoundarySlabSet(int j_decay, multi1d<int> ti,multi1d<int> th,int N) : decay_dir(j_decay),t0(ti),thick(th),Nt(N)
  {
    sets.make(FixedBoundaryFunc(j_decay,t0,thick,N));
  }


 

}  // end namespace Chroma
