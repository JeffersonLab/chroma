// -*- C++ -*-
// $Id: invbicgstab.h,v 1.1 2004-05-19 00:21:23 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __invbicgstab__
#define __invbicgstab__

#include "linearop.h"

template<typename T>
void InvBiCGStab(const LinearOperator<T>& A,
		 const T& chi,
		 T& psi,
		 const Real& RsdCG, 
		 int MaxCG, 
		 int& n_count);

#endif
