// -*- C++ -*-
// $Id: enum_loops_s.h,v 2.1 2006-01-31 14:31:34 mcneile Exp $
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
    GAUSSIAN,
    T_DILUTE_GAUSS,
    C_DILUTE_GAUSS,
    P_DILUTE_GAUSS,
    CT_DILUTE_GAUSS,
    CP_DILUTE_GAUSS,
    PT_DILUTE_GAUSS,
    MOD_T_DILUTE_GAUSS,
    CORNER_DILUTE_GAUSS,
    COR_DBL_T_DILUTE_GAUSS,
    COR_MOD_DBL_T_DILUTE_GAUSS,
    C_MOD_T_DILUTE_GAUSS 
  };
  typedef   VolSrc  VolSrc_type ;




}  // end namespace Chroma

#endif
