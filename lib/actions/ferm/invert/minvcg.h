// -*- C++ -*-
// $Id: minvcg.h,v 1.1 2003-04-08 21:36:13 edwards Exp $

#ifndef MINVCG_INCLUDE
#define MINVCG_INCLUDE

void MInvCG(const LinearOperator& A, const LatticeFermion& chi, LatticeFermion& psi,
	    const Double& chi_norm, 
	    const multi1d<Real>& mass, int Nmass, int isz, 
	    const multi1d<Real>& RsdCG, int& n_count,
	    const multi1d<LatticeFermion>& EigVec, int NEig);

#endif
