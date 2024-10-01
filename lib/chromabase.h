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

#include <map>
#include <exception>
#include <typeinfo>
#include <cassert>

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

#if defined (QDP_IS_QDPJIT2)
template<>
struct PropTypeTraits<LatticeFermion>
{
  typedef LatticePropagator   Type_t;
};
#endif

template<>
struct PropTypeTraits<LatticeStaggeredFermion>
{
  typedef LatticeStaggeredPropagator   Type_t;
};


// Extremely basic types
enum PlusMinus {PLUS = 1, MINUS = -1};


  struct __chroma_constant
  {
    Real twopi;
    Real fuzz;
  };
  const __chroma_constant& constant();
  void constant_destroy();
  

// Hooks for various things
#if defined(QDP_DEBUG_MEMORY)
#define START_CODE() QDP::Allocator::theQDPAllocator::Instance().pushFunc(__func__, __LINE__)
#define END_CODE()   QDP::Allocator::theQDPAllocator::Instance().popFunc()

#else
#define START_CODE() QDP_PUSH_PROFILE(QDP::getProfileLevel())
#define END_CODE()   QDP_POP_PROFILE()

#endif
}

#endif
