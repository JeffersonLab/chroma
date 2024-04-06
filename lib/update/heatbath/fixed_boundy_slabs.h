// -*- C++ -*-
/*! \file
 *  \brief Convenience for building time-slab subsets
 */

#ifndef __fixed_bc_slabs_set_h__
#define __fixed_bc_slabs_set_h__

#include "chromabase.h"

namespace Chroma 
{

  //! Builds time slab subsets
  /*!
   * \ingroup ft
   */
  class FixedBoundarySlabSet
  {
  public:
    //! Constructor about origin
    FixedBoundarySlabSet(int decay_dir,multi1d<int> t0, multi1d<int> thick,int N);
    
    //! The set to be used in sumMulti
    const Set& getSet() const {return sets;}

    //! Number of subsets - There are 1 or 2 subsets 
    int numSubsets() const {return sets.numSubsets();}

    //! Decay direction
    int getDir() const {return decay_dir;}
    
    multi1d<int> getT0() const {return t0 ;}

    multi1d<int> getThick() const {return thick ;}

    int getNt() const {return Nt ;}
  private:
    FixedBoundarySlabSet() {} // hide default constructor

    int  decay_dir;
    multi1d<int>  t0 ;
    multi1d<int>  thick ;
    int  Nt ;
    Set  sets;
  };
  
}  // end namespace Chroma

#endif
