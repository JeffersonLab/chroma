// -*- C++ -*-
// $Id: enum_loops_s.h,v 1.1 2005-08-27 11:34:08 mcneile Exp $
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
