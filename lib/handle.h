// -*- C++ -*-
/*! @file
 * @brief  Class for counted reference semantics
 *
 * Holds and object, and deletes it when the last Handle to it
 * is destroyed
 *
 * Code from  * "The C++ Standard Library - A Tutorial and Reference"
 * by Nicolai M. Josuttis, Addison-Wesley, 1999.
 *
 * An almost identical version is in Stroustrup, "C++ Programming Language",
 * 3rd ed., Section 25.7. The names used there are used for this class.
 *
 */

#ifndef __handle_h__
#define __handle_h__

#include "chromabase.h"
#include <iostream>

using namespace QDP;

namespace Chroma
{
  //! Class for counted reference semantics
  /*!
   * Holds and object, and deletes it when the last Handle to it
   * is destroyed
   */


  template<typename T>
  class Handle
  {
  public:

    Handle(Handle&& other) noexcept : ptr(other.ptr), count(other.count)
    {
      other.ptr = nullptr;
      other.count = nullptr; // Or a new int(0) depending on your design
    }

    Handle& operator=(Handle&& other) noexcept
    {
      if (this != &other)
	{
	  dispose();

	  ptr = other.ptr;
	  count = other.count;

	  other.ptr = nullptr;
	  other.count = nullptr;
	}
      return *this;
    }
    

    Handle(T* p = nullptr) : ptr(p) {
      if (ptr) {
	count = new int(1);
      } else {
	count = new int(0); 
      }
    }

    Handle(const Handle& p) : ptr(p.ptr), count(p.count)
    {
      if (ptr) {
	++*count;
      }
    }

    ~Handle() { dispose(); }

    Handle& operator=(const Handle& p)
    {
      if (this != &p)
	{
	  dispose();
	  ptr = p.ptr;
	  count = p.count;
	  if (ptr) {
	    ++*count;
	  }
	}
      return *this;
    }


    template<typename Q>
    Handle<Q> cast() const
    {
      Handle<Q> q = this->tryCast<Q>();
      if (!q) {
	QDPIO::cerr << "Dynamic cast failed in Handle::cast()" << std::endl;
	QDPIO::cerr << "You are trying to cast to a class you cannot cast to" << std::endl;
	QDP_abort(1);
      }
      return q;
    }

    
    template<typename Q>
    Handle<Q> tryCast() const
    {
      if (Q* new_ptr = dynamic_cast<Q*>(ptr)) {
	return Handle<Q>(new_ptr, count);
      }
      return Handle<Q>();
    }


    explicit operator bool() const { return ptr != nullptr; }


    template<typename Q> friend class Handle;


    T& operator*() const { assert(ptr); return *ptr; }
    T* operator->() const { assert(ptr); return ptr; }
    T* get() const { return ptr; }

  private:
    Handle(T* p, int* c) : ptr(p), count(c)
    {
      if (ptr) {
	++*count;
      }
    }

    void dispose()
    {
      if (ptr && --*count == 0)
	{
	  delete count;
	  delete ptr;
	}
    }

  private:
    T* ptr;   // Pointer to the value
    int* count; // Shared number of owners
  };



}// namespace


#endif
