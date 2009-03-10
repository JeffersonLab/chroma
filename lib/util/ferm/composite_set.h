// -*- C++ -*-
// $Id: composite_set.h,v 1.2 2009-03-10 21:55:16 kostas Exp $
/*! \file
 * \brief Key for propagator colorvector sources
 */

#ifndef __FUNC_COMPOSITE_SET_h__
#define __FUNC_COMPOSITE_SET_h__

#include "chromabase.h"

namespace Chroma
{

  //! Function object used for constructing a composite  set
  class CompositeSetFunc : public SetFunc
  {
  public:
    CompositeSetFunc(const SetFunc& s1, const SetFunc& s2 ):set1(s1),set2(s2){
    }

    int operator() (const multi1d<int>& coordinate) const{
      return set1(coordinate) + set1.numSubsets()*set2(coordinate) ;
    }

    int numSubsets() const{
      return set1.numSubsets()*set2.numSubsets() ;
    }

    int coloring(int s1, int s2){
      return s1 + set1.numSubsets()*s2 ;
    }

  private:
    CompositeSetFunc() {}  // hide default constructor
    
    SetFunc set1 ;
    SetFunc set2 ;

  };

} // namespace Chroma

#endif
