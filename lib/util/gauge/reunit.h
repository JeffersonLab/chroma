// -*- C++ -*-
// $Id: reunit.h,v 1.1 2003-02-15 05:54:26 edwards Exp $

#ifndef REUNIT_INCLUDE
#define REUNIT_INCLUDE

enum Reunitarize {REUNITARIZE, REUNITARIZE_ERROR, REUNITARIZE_LABEL};

void reunit(LatticeColorMatrix& xa);
void reunit(LatticeColorMatrix& xa, LatticeBoolean& bad, int& numbad, enum Reunitarize ruflag);

#endif
