// -*- C++ -*-
// $Id: z2_src.h,v 1.4 2005-01-14 18:42:37 edwards Exp $
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
