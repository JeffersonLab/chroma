// -*- C++ -*-
// $Id: handle.h,v 3.4 2007-08-12 20:35:11 edwards Exp $
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
using namespace std;

namespace Chroma
{
  //! Class for counted reference semantics
  /*!
   * Holds and object, and deletes it when the last Handle to it
   * is destroyed
   */
  template <class T>
  class Handle
  {
  public:
    //! Initialize pointer with existing pointer
    /*! Requires that the pointer p is a return value of new */
    Handle(T* p=0) : ptr(p), count(new int(1)) {}

    //! Copy pointer (one more owner)
    Handle(const Handle& p) : ptr(p.ptr), count(p.count) 
      {++*count;}

    //! Destructor (delete value if this was the last owner)
    ~Handle() {dispose();}

    //! Assignment (unshare old and share new value)
    Handle& operator=(const Handle& p) 
      {
	if (ptr != p.ptr)
	{
	  dispose();
	  ptr = p.ptr;
	  count = p.count;
	  ++*count;
	}
	return *this;
      }

    //! RGE's addition. A cast function to morph the actual type
    template<typename Q>
    Handle<Q> cast()
      {
	Handle<Q> q;
	q.ptr = dynamic_cast<Q*>(ptr);
	if( q.ptr == 0x0 ) { 
	  QDPIO::cerr << "Dynamic cast failed in Handle::cast()" <<endl;
	  QDPIO::cerr << "You are trying to cast to a class you cannot cast to" << endl;
	  QDP_abort(1);
	}
	delete q.count;

	q.count = count;
	++*count;
	return q;
      }

    //! The cast function requires all Handles<Q> to be friends of Handle<T>
    template<typename Q> friend class Handle;

    //! Access the value to which the pointer refers
    T& operator*() const {return *ptr;}
    T* operator->() const {return ptr;}

  private:
    void dispose() 
      {
	if (--*count == 0) 
	{
	  delete count;
	  delete ptr;
	}
      }

  private:
    T* ptr;        // pointer to the value
    int* count;    // shared number of owners
  };

}


#endif
