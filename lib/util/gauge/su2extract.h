// -*- C++ -*-
// $Id: su2extract.h,v 1.1 2003-03-28 05:34:39 edwards Exp $

#ifndef __su2extract__
#define __su2extract__

multi1d<LatticeReal> 
su2Extract(const LatticeColorMatrix& source, 
	   int su2_index, 
	   const Subset& s);

#endif
