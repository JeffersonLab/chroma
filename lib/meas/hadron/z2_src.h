// -*- C++ -*-
// $Id: z2_src.h,v 1.3 2004-11-20 19:25:10 mcneile Exp $
/*! \file
 *  \brief Volume source of Z2 noise
 */


//
//
//

#ifndef  Z2_SRC_INC
#define  Z2_SRC_INC 

void z2_src(LatticeFermion& a) ;
void z2_src(LatticeStaggeredFermion& a) ; 

void z2_src(LatticeFermion& a, int slice, int mu) ;

#endif
