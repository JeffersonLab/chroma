// -*- C++ -*-
// $Id: chroma_scalarsite_defs.h,v 1.2 2004-03-29 21:31:15 edwards Exp $

/*! \file
 * \brief Additional QDP type definitions
 */

#include "qdp.h"

QDP_BEGIN_NAMESPACE(QDP);

/*! \addtogroup defs Type definitions
 *
 * User constructed types made from QDP type compositional nesting.
 * The layout is suitable for a scalar-like implementation. Namely,
 * sites are the slowest varying index.
 *
 * @{
 */

const int Ls = 4;   // HACK - FOR now fixed DW index here

// Aliases for a scalar architecture

typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, Nc>, Ns> > LatticeDWFermion;
typedef OLattice< PSpinVector< PColorVector< RComplex<REAL>, Nc>, Ns>>1> > LatticeHalfDWFermion;

typedef OScalar< PSpinVector< PColorVector< RComplex<REAL>, Nc>, Ns> > DWFermion;


/*! @} */   // end of group defs

QDP_END_NAMESPACE();

