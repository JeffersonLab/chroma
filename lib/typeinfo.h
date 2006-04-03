// -*- C++ -*-
// $Id: typeinfo.h,v 3.0 2006-04-03 04:58:44 edwards Exp $
/*! @file
 * @brief Type info support
 */

#ifndef __typeinfo_h__
#define __typeinfo_h__

#include <typeinfo>
#include <cassert>
#include "typelist.h"

namespace Chroma
{
////////////////////////////////////////////////////////////////////////////////
// class TypeInfo
// Purpose: offer a first-class, comparable wrapper over std::type_info
////////////////////////////////////////////////////////////////////////////////

  class TypeInfo
  {
  public:
    // Constructors
    TypeInfo(); // needed for containers
    TypeInfo(const std::type_info&); // non-explicit

    // Access for the wrapped std::type_info
    const std::type_info& Get() const;
    // Compatibility functions
    bool before(const TypeInfo& rhs) const;
    const char* name() const;

  private:
    const std::type_info* pInfo_;
  };
    
// Implementation
    
  inline TypeInfo::TypeInfo()
  {
    class Nil {};
    pInfo_ = &typeid(Nil);
    assert(pInfo_);
  }
    
  inline TypeInfo::TypeInfo(const std::type_info& ti)
    : pInfo_(&ti)
  { assert(pInfo_); }
    
  inline bool TypeInfo::before(const TypeInfo& rhs) const
  {
    assert(pInfo_);
    return pInfo_->before(*rhs.pInfo_);
  }

  inline const std::type_info& TypeInfo::Get() const
  {
    assert(pInfo_);
    return *pInfo_;
  }
    
  inline const char* TypeInfo::name() const
  {
    assert(pInfo_);
    return pInfo_->name();
  }

// Comparison operators
    
  inline bool operator==(const TypeInfo& lhs, const TypeInfo& rhs)
  { return lhs.Get() == rhs.Get(); }

  inline bool operator<(const TypeInfo& lhs, const TypeInfo& rhs)
  { return lhs.before(rhs); }

  inline bool operator!=(const TypeInfo& lhs, const TypeInfo& rhs)
  { return !(lhs == rhs); }    
    
  inline bool operator>(const TypeInfo& lhs, const TypeInfo& rhs)
  { return rhs < lhs; }
    
  inline bool operator<=(const TypeInfo& lhs, const TypeInfo& rhs)
  { return !(lhs > rhs); }
     
  inline bool operator>=(const TypeInfo& lhs, const TypeInfo& rhs)
  { return !(lhs < rhs); }
}

#endif

