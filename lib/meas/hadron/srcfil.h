// -*- C++ -*-
// $Id: srcfil.h,v 1.1 2003-02-15 05:54:26 edwards Exp $

#ifndef SRCFIL_INCLUDE
#define SRCFIL_INCLUDE

void srcfil(LatticeFermion& a, multi1d<int>& coord, int color_index, int spin_index);

void srcfil(LatticeColorVector& a, multi1d<int>& coord, int color_index);

#endif
