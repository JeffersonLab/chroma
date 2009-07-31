#ifndef POINTER_H
#define POINTER_H

#include <iostream>
#include <list>
#include "qdp.h"
#include "chromabase.h"

namespace Chroma {
  namespace LaphEnv {

// ******************************************************************
// *                                                                *
// *  Pointers are a powerful feature of C++, but they are          *
// *  prone to errors and debugging difficulties.  The purpose      *
// *  of the class "Pointer" and "constPointer" below is to         *
// *  provide a safer way of pointing to certain objects.           *
// *  In particular, if the object pointed to goes out of scope,    *
// *  the "Pointer" to it knows and can be tested for nullness.     *
// *                                                                *
// *  To provide "Pointer" objects to a class "T", you must do      *
// *  the following:                                                *
// *                                                                *
// *  (1) In the "private" part of class "T", include the lines     *
// *                                                                *
// *     mutable list<Pointer<T>* > ptrlist;   (replace T by        *
// *     friend class Pointer<T>;               the class name)     *
// *                                                                *
// *  (2) In the destructor of class "T", include a statement       *
// *                                                                *
// *     Pointer<T>::nullify_pointers(ptrlist);                     *
// *                                                                *
// *  That's it!  Then you can create Pointer<T> objects and        *
// *  use such objects just like ordinary pointers.                 *
// *                                                                *
// *  Usage:                                                        *
// *                                                                *
// *    class Foo { ... }  // defined as above                      *
// *                                                                *
// *    Foo anobj(...);                                             *
// *    Pointer<Foo> fptr(&anobj);                                  *
// *    fptr->foo_member(...);                                      *
// *                                                                *
// *    bool b = fptr.isNull();  // true if fptr points to nothing  *
// *    b = fptr.isNotNull();    // true if fptr points to valid    *
// *                                                                *
// *                                                                *
// *  To provide "constPointer" objects to a class "T", you must    *
// *  do the following:                                             *
// *                                                                *
// *  (1) In the "private" part of class "T", include the lines     *
// *                                                                *
// *     mutable list<constPointer<T>* > constptrlist; (replace T by*
// *     friend class constPointer<T>;               the class name)*
// *                                                                *
// *  (2) In the destructor of class "T", include a statement       *
// *                                                                *
// *     constPointer<T>::nullify_pointers(constptrlist);           *
// *                                                                *
// *                                                                *
// *                                                                *
// *  Implementation notes:  The "ptrlist" of class "T" stores all  *
// *  Pointer<T> objects that point to it.  When the "T" object     *
// *  goes out of scope, the "nullify_pointers" subroutine first    *
// *  nullifies all Pointer<T> objects pointing to the object first.*
// *  Creating and deleting Pointer<T> objects either adds to or    *
// *  removes these Pointer<T> objects from the "ptrlist".          *
// *                                                                *
// *  Reminder:  const Fred* p; // p can be changed, but not Fred   *
// *             Fred* const p; // p constant, but Fred changeable  *
// *                                                                *
// ******************************************************************


template <typename T>
class Pointer
{
    T* ptr;

  public:

    Pointer() : ptr(0) {}

    Pointer(T *const anobj);

    Pointer(const Pointer<T>& aptr);
    
    ~Pointer();

    Pointer& operator=(T *const anobj);

    Pointer& operator=(const Pointer<T>& aptr);

    T* operator->() const
    {
     if (ptr==0){
        QDPIO::cerr << "attempt to access null Pointer"<<std::endl;
        QDP_abort(1);}
     return ptr;
    }

    T& operator*() const
    {
     if (ptr==0){
        QDPIO::cerr << "attempt to access null Pointer"<<std::endl;
        QDP_abort(1);}
     return *ptr;
    }

    bool isNull() const { return (ptr==0); }

    bool isNotNull() const { return (ptr!=0); }

  private:

   void add_pointer_to_object();
   void remove_pointer_from_object();
   static void nullify_pointers(std::list<Pointer<T>* >& ptrlist);

   friend T::~T();

};


template <typename T>
Pointer<T>::Pointer(T *const anobj)
{
 ptr=anobj;
 add_pointer_to_object();
}

template <typename T>
Pointer<T>::Pointer(const Pointer<T>& aptr)
{
 ptr=aptr.ptr;
 add_pointer_to_object();
}
    

template <typename T>
Pointer<T>& Pointer<T>::operator=(T *const anobj)
{
 if (ptr!=0) remove_pointer_from_object();
 ptr=anobj;
 add_pointer_to_object();
}


template <typename T>
Pointer<T>& Pointer<T>::operator=(const Pointer<T>& aptr)
{
 if (ptr!=0) remove_pointer_from_object();
 ptr=aptr.ptr;
 add_pointer_to_object();
}


template <typename T>
Pointer<T>::~Pointer()
{
 if (ptr!=0) remove_pointer_from_object();
}


template <typename T>
void Pointer<T>::add_pointer_to_object()
{
 ptr->ptrlist.push_back(this);
}

template <typename T>
void Pointer<T>::remove_pointer_from_object()
{
 ptr->ptrlist.remove(this);
}

template <typename T>
void Pointer<T>::nullify_pointers(std::list<Pointer<T>* >& ptrlist)
{
 typename std::list< Pointer<T>* >::iterator it;
 for (it=ptrlist.begin();it!=ptrlist.end();it++){
    (*it)->ptr=0;}
 ptrlist.clear(); 
}


 // *****************************************************************

template <typename T>
class constPointer
{
    const T* ptr;

  public:

    constPointer() : ptr(0) {}

    constPointer(const T *const anobj);

    constPointer(const constPointer<T>& aptr);
    
    ~constPointer();

    constPointer& operator=(const T *const anobj);

    constPointer& operator=(const constPointer<T>& aptr);

    const T* operator->() const
    {
     if (ptr==0){
        QDPIO::cerr << "attempt to access null constPointer"<<std::endl;
        QDP_abort(1);}
     return ptr;
    }

    const T& operator*() const
    {
     if (ptr==0){
        QDPIO::cerr << "attempt to access null constPointer"<<std::endl;
        QDP_abort(1);}
     return *ptr;
    }

    bool isNull() const { return (ptr==0); }

    bool isNotNull() const { return (ptr!=0); }

  private:

   void add_pointer_to_object();
   void remove_pointer_from_object();
   static void nullify_pointers(std::list<constPointer<T>* >& constptrlist);

   friend T::~T();

};


template <typename T>
constPointer<T>::constPointer(const T *const anobj)
{
 ptr=anobj;
 add_pointer_to_object();
}

template <typename T>
constPointer<T>::constPointer(const constPointer<T>& aptr)
{
 ptr=aptr.ptr;
 add_pointer_to_object();
}
    

template <typename T>
constPointer<T>& constPointer<T>::operator=(const T *const anobj)
{
 if (ptr!=0) remove_pointer_from_object();
 ptr=anobj;
 add_pointer_to_object();
}


template <typename T>
constPointer<T>& constPointer<T>::operator=(const constPointer<T>& aptr)
{
 if (ptr!=0) remove_pointer_from_object();
 ptr=aptr.ptr;
 add_pointer_to_object();
}


template <typename T>
constPointer<T>::~constPointer()
{
 if (ptr!=0) remove_pointer_from_object();
}


template <typename T>
void constPointer<T>::add_pointer_to_object()
{
 ptr->constptrlist.push_back(this);
}

template <typename T>
void constPointer<T>::remove_pointer_from_object()
{
 ptr->constptrlist.remove(this);
}

template <typename T>
void constPointer<T>::nullify_pointers(std::list<constPointer<T>* >& ptrlist)
{
 typename std::list< constPointer<T>* >::iterator it;
 for (it=ptrlist.begin();it!=ptrlist.end();it++){
    (*it)->ptr=0;}
 ptrlist.clear(); 
}


// *************************************************************
  }
}
#endif
