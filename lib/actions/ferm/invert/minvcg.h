// -*- C++ -*-
// $Id: minvcg.h,v 1.3 2003-10-20 20:31:50 edwards Exp $

#ifndef MINVCG_INCLUDE
#define MINVCG_INCLUDE

template<typename T>
void MInvCG(const LinearOperator<T>& A, 
	    const T& chi, 
	    T& psi,
	    const Double& chi_norm, 
	    const multi1d<Real>& mass, int Nmass, int isz, 
	    const multi1d<Real>& RsdCG, int& n_count,
	    const multi1d<T>& EigVec, int NEig);

#endif
