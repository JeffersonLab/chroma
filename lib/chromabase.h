// $Id: chromabase.h,v 1.3 2004-07-27 06:00:07 edwards Exp $
//
// Absolute basic stuff to use chroma
//
/*! \file
 * \brief Primary include file for CHROMA library code
 *
 * This is the absolute basic stuff to use Chroma in say
 * library codes.
 */

#ifndef CHROMABASE_INCLUDE
#define CHROMABASE_INCLUDE

#include "qdp.h"
#include "chroma_config.h"

namespace QDP {

// Trait classes to undo DW index
template<class T>
struct BaseType {};

template<>
struct BaseType<LatticeFermion>
{
  typedef LatticeFermion   Type_t;
};
}

using namespace QDP;

// Extremely basic types
enum PlusMinus {PLUS = 1, MINUS = -1};


// Useful constants
#if BASE_PRECISION == 32
const Real fuzz = 1.0e-5;
#elif BASE_PRECISION == 64
const Real fuzz = 1.0e-10;
#endif

const Real twopi = 6.283185307179586476925286;

#define TO_REAL(a) float(a)


// Hooks for various things
// NOTE: can get rid of unused arg "a" !!
#define START_CODE(a) QDP_PUSH_PROFILE(QDP::getProfileLevel())
#define END_CODE(a)   QDP_POP_PROFILE()


#endif
