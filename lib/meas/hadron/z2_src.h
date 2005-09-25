// -*- C++ -*-
// $Id: z2_src.h,v 2.0 2005-09-25 21:04:36 edwards Exp $
/*! \file
 *  \brief Volume source of Z2 noise
 */


//
//
//

#ifndef  Z2_SRC_INC
#define  Z2_SRC_INC 

namespace Chroma {

void z2_src(LatticeFermion& a) ;
void z2_src(LatticeStaggeredFermion& a) ; 

void z2_src(LatticeFermion& a, int slice, int mu) ;

}  // end namespace Chroma

#endif
