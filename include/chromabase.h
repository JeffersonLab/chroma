// $Id: chromabase.h,v 1.7 2003-10-20 20:23:48 edwards Exp $
//
// Absolute basic stuff to use chroma

#ifndef CHROMABASE_INCLUDE
#define CHROMABASE_INCLUDE

#include "qdp.h"
#include "qdp_util.h"

#if defined(ARCH_SCALAR) || defined(ARCH_PARSCALAR)
#include "chroma_scalarsite_defs.h"

#elif defined(ARCH_SCALARVEC) || defined(ARCH_PARSCALARVEC)
// #include "chroma_scalarvecsite_defs.h"

#else
#error "Unknown architecture ARCH"
#endif



using namespace QDP;

const float fuzz = 1.0e-5;
const float twopi = 6.283185307179586476925286;

#define START_CODE(a)
#define END_CODE(a)

#define TO_REAL(a) float(a)

#endif
