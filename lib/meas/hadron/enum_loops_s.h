// -*- C++ -*-
// $Id: enum_loops_s.h,v 2.0 2005-09-25 21:04:35 edwards Exp $
/*! \file
 *  \brief Enums for the different types of stochastic source.
 */        


//
//
//

#ifndef  ENUM_LOOPS_S_INC
#define  ENUM_LOOPS_S_INC

namespace Chroma {

  enum VolSrc {
    Z2NOISE ,
    GAUSSIAN
  };
  typedef   VolSrc  VolSrc_type ;




}  // end namespace Chroma

#endif
