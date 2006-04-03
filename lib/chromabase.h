// $Id: chromabase.h,v 3.0 2006-04-03 04:58:43 edwards Exp $
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
// #include "chroma_config.h"   // turn off by default using config file

using namespace QDP;

namespace Chroma {

// Trait classes 
template<class T>
struct PropTypeTraits {};

template<>
struct PropTypeTraits<LatticeDiracFermion>
{
  typedef LatticeDiracPropagator   Type_t;
};


template<>
struct PropTypeTraits<LatticeStaggeredFermion>
{
  typedef LatticeStaggeredPropagator   Type_t;
};


// Extremely basic types
enum PlusMinus {PLUS = 1, MINUS = -1};


// Useful constants
#if BASE_PRECISION == 32
const Real fuzz = 1.0e-5;
#elif BASE_PRECISION == 64
const Real fuzz = 1.0e-10;
#endif

const Real twopi = 6.283185307179586476925286;


// Hooks for various things
#if defined(QDP_DEBUG_MEMORY)
#define START_CODE() QDP::Allocator::theQDPAllocator::Instance().pushFunc(__func__, __LINE__)
#define END_CODE()   QDP::Allocator::theQDPAllocator::Instance().popFunc()

#else
#define START_CODE() QDP_PUSH_PROFILE(QDP::getProfileLevel())
#define END_CODE()   QDP_POP_PROFILE()

#endif
};

#endif
