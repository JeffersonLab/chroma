// -*- C++ -*-
/*! @file
 * @brief Type manipulator support
 */

#ifndef __typemanip_h__
#define __typemanip_h__

namespace Chroma
{

////////////////////////////////////////////////////////////////////////////////
// class template Int2Type
// Converts each integral constant into a unique type
// Invocation: Int2Type<v> where v is a compile-time constant integral
// Defines 'value', an enum that evaluates to v
////////////////////////////////////////////////////////////////////////////////

  template <int v>
  struct Int2Type
  {
    enum { value = v };
  };
    
////////////////////////////////////////////////////////////////////////////////
// class template Type2Type
// Converts each type into a unique, insipid type
// Invocation Type2Type<T> where T is a type
// Defines the type OriginalType which std::maps back to T
////////////////////////////////////////////////////////////////////////////////

  template <typename T>
  struct Type2Type
  {
    typedef T OriginalType;
  };
    
////////////////////////////////////////////////////////////////////////////////
// class template Select
// Selects one of two types based upon a boolean constant
// Invocation: Select<flag, T, U>::Result
// where:
// flag is a compile-time boolean constant
// T and U are types
// Result evaluates to T if flag is true, and to U otherwise.
////////////////////////////////////////////////////////////////////////////////

  template <bool flag, typename T, typename U>
  struct Select
  {
    typedef T Result;
  };
  template <typename T, typename U>
  struct Select<false, T, U>
  {
    typedef U Result;
  };
    
////////////////////////////////////////////////////////////////////////////////
// Helper types Small and Big - guarantee that sizeof(Small) < sizeof(Big)
////////////////////////////////////////////////////////////////////////////////

  namespace Private
  {
    template <class T, class U>
    struct ConversionHelper
    {
      typedef char Small;
      struct Big { char dummy[2]; };
      static Big   Test(...);
      static Small Test(U);
      static T MakeT();
    };
  }

////////////////////////////////////////////////////////////////////////////////
// class template Conversion
// Figures out the conversion relationships between two types
// Invocations (T and U are types):
// a) Conversion<T, U>::exists
// returns (at compile time) true if there is an implicit conversion from T
// to U (example: Derived to Base)
// b) Conversion<T, U>::exists2Way
// returns (at compile time) true if there are both conversions from T
// to U and from U to T (example: int to char and back)
// c) Conversion<T, U>::sameType
// returns (at compile time) true if T and U represent the same type
//
// Caveat: might not work if T and U are in a private inheritance hierarchy.
////////////////////////////////////////////////////////////////////////////////

  template <class T, class U>
  struct Conversion
  {
    typedef Private::ConversionHelper<T, U> H;
#ifndef __MWERKS__
    enum { exists = sizeof(typename H::Small) == sizeof(H::Test(H::MakeT())) };
#else
    enum { exists = false };
#endif
    enum { exists2Way = exists && Conversion<U, T>::exists };
    enum { sameType = false };
  };
    
  template <class T>
  struct Conversion<T, T>    
  {
    enum { exists = 1, exists2Way = 1,sameType = 1 };
  };
    
  template <class T>
  struct Conversion<void, T>    
  {
    enum { exists = 1, exists2Way = 0,sameType = 0 };
  };
    
  template <class T>
  struct Conversion<T, void>    
  {
    enum { exists = 1, exists2Way = 0,sameType = 0 };
  };
    
  template <>
  class Conversion<void, void>    
  {
  public:
    enum { exists = 1, exists2Way = 1,sameType = 1 };
  };    
}

////////////////////////////////////////////////////////////////////////////////
// macro SUPERSUBCLASS
// Invocation: SUPERSUBCLASS(B, D) where B and D are types. 
// Returns true if B is a public base of D, or if B and D are aliases of the 
// same type.
//
// Caveat: might not work if T and U are in a private inheritance hierarchy.
////////////////////////////////////////////////////////////////////////////////

#define SUPERSUBCLASS(T, U) \
    (::Chroma::Conversion<const U*, const T*>::exists && \
    !::Chroma::Conversion<const T*, const void*>::sameType)

////////////////////////////////////////////////////////////////////////////////
// macro SUPERSUBCLASS
// Invocation: SUPERSUBCLASS(B, D) where B and D are types. 
// Returns true if B is a public base of D.
//
// Caveat: might not work if T and U are in a private inheritance hierarchy.
////////////////////////////////////////////////////////////////////////////////

#define SUPERSUBCLASS_STRICT(T, U) \
    (SUPERSUBCLASS(T, U) && \
    !::Chroma::Conversion<const T, const U>::sameType)

#endif
