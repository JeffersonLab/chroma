// -*- C++ -*-
// $Id: sunfill.h,v 1.3 2003-08-09 03:56:11 edwards Exp $

#ifndef __sunfill_h__
#define __sunfill_h__

LatticeColorMatrix
sunFill(const multi1d<LatticeReal>& r,
	int su2_index,
	const UnorderedSubset& s);

LatticeColorMatrix
sunFill(const multi1d<LatticeReal>& r,
	int su2_index,
	const OrderedSubset& s);

#endif
