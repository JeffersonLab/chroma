// -*- C++ -*-
// $Id: sunfill.h,v 1.2 2003-04-02 22:24:58 edwards Exp $

#ifndef __sunfill_h__
#define __sunfill_h__

LatticeColorMatrix
sunFill(const multi1d<LatticeReal>& r,
	int su2_index,
	const Subset& s);

#endif
