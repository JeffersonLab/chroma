// -*- C++ -*-
// $Id: invbicgstab_array.h,v 1.1 2004-05-19 00:21:23 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __invbicgstab_array__
#define __invbicgstab_array__

#include "linearop.h"

template<typename T>
void InvBiCGStab(const LinearOperator< multi1d<T> >& A,
		 const multi1d<T>& chi,
		 multi1d<T>& psi,
		 const Real& RsdCG, 
		 int MaxCG, 
		 int& n_count);

#endif
