// -*- C++ -*-
// $Id: chroma_scalarsite_dwdefs.h,v 1.2 2003-10-22 03:22:09 edwards Exp $

/*! \file
 * \brief Additional QDP type definitions
 */

#include "qdp.h"

#define CHROMA_USE_DWF

QDP_BEGIN_NAMESPACE(QDP);

/*! \addtogroup defs Type definitions
 *
 * User constructed types made from QDP type compositional nesting.
 * The layout is suitable for a scalar-like implementation. Namely,
 * sites are the slowest varying index.
 *
 * @{
 */

const int Ls = 8;   // HACK - FOR now fixed DW index here

// Aliases for a scalar architecture

// New guy beyond qdp_defs.h
typedef OLattice< PDWVector< PSpinVector< PColorVector< RComplex<REAL>, Nc>, Ns>, Ls> > LatticeDWFermion;
typedef OLattice< PDWVector< PSpinVector< PColorVector< RComplex<REAL>, Nc>, Ns>>1>, Ls> > LatticeHalfDWFermion;

typedef OScalar< PDWVector< PSpinVector< PColorVector< RComplex<REAL>, Nc>, Ns>, Ls> > DWFermion;


/*! @} */   // end of group defs

QDP_END_NAMESPACE();

