// $Id: chromabase.h,v 1.9 2003-11-09 22:36:13 edwards Exp $
//
// Absolute basic stuff to use chroma

#ifndef CHROMABASE_INCLUDE
#define CHROMABASE_INCLUDE

#include "qdp.h"

#if defined(ARCH_SCALAR) || defined(ARCH_PARSCALAR)
#include "chroma_scalarsite_dwdefs.h"

#elif defined(ARCH_SCALARVEC) || defined(ARCH_PARSCALARVEC)
// #include "chroma_scalarvecsite_defs.h"

#else
#error "Unknown architecture ARCH"
#endif

using namespace QDP;

// Extremely basic types
enum PlusMinus {PLUS = 1, MINUS = -1};


// Useful constants
const float fuzz = 1.0e-5;
const float twopi = 6.283185307179586476925286;

#define TO_REAL(a) float(a)


// Hooks for various things
#define START_CODE(a)
#define END_CODE(a)


#endif
