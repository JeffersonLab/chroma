// -*- C++ -*-
// $Id: su2extract.h,v 1.2 2003-08-09 03:56:11 edwards Exp $

#ifndef __su2extract__
#define __su2extract__

multi1d<LatticeReal> 
su2Extract(const LatticeColorMatrix& source, 
	   int su2_index, 
	   const UnorderedSubset& s);

multi1d<LatticeReal> 
su2Extract(const LatticeColorMatrix& source, 
	   int su2_index, 
	   const OrderedSubset& s);

#endif
