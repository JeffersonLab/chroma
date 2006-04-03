// -*- C++ -*-
// $Id: objfunctor.h,v 3.0 2006-04-03 04:58:44 edwards Exp $
/*! @file
 * @brief Generic functor class
 */

#ifndef __objfunctor_h__
#define __objfunctor_h__

#include "typelist.h"
#include "typetraits.h"
#include <typeinfo>
#include <memory>

namespace Chroma
{
////////////////////////////////////////////////////////////////////////////////
// class template FunctorImpl (internal)
////////////////////////////////////////////////////////////////////////////////

  namespace Private
  {
    template <typename R>
    struct FunctorImplBase 
    {
      typedef R ResultType;
            
      typedef EmptyType Parm1;
      typedef EmptyType Parm2;
      typedef EmptyType Parm3;
      typedef EmptyType Parm4;
      typedef EmptyType Parm5;

      virtual FunctorImplBase* DoClone() const = 0;
      template <class U>
      static U* Clone(U* pObj)
	{
	  if (!pObj) return 0;
	  U* pClone = static_cast<U*>(pObj->DoClone());
	  assert(typeid(*pClone) == typeid(*pObj));
	  return pClone;
	}
    };
  }
    
////////////////////////////////////////////////////////////////////////////////
// macro DEFINE_CLONE_FUNCTORIMPL
// Implements the DoClone function for a functor implementation
////////////////////////////////////////////////////////////////////////////////

#define DEFINE_CLONE_FUNCTORIMPL(Cls) \
    virtual Cls* DoClone() const { return new Cls(*this); }

////////////////////////////////////////////////////////////////////////////////
// class template FunctorImpl
// The base class for a hierarchy of functors. The FunctorImpl class is not used
//     directly; rather, the Functor class manages and forwards to a pointer to
//     FunctorImpl
// You may want to derive your own functors from FunctorImpl.
// Specializations of FunctorImpl for up to 5 parameters follow
////////////////////////////////////////////////////////////////////////////////

  template <typename R, class TList>
  class FunctorImpl;

////////////////////////////////////////////////////////////////////////////////
// class template FunctorImpl
// Specialization for 0 (zero) parameters
////////////////////////////////////////////////////////////////////////////////

  template <typename R>
  class FunctorImpl<R, NullType>
    : public Private::FunctorImplBase<R>
  {
  public:
    typedef R ResultType;
    virtual R operator()() = 0;
  };

////////////////////////////////////////////////////////////////////////////////
// class template FunctorImpl
// Specialization for 1 parameter
////////////////////////////////////////////////////////////////////////////////

  template <typename R, typename P1>
  class FunctorImpl<R, TYPELIST_1(P1)>
    : public Private::FunctorImplBase<R>
  {
  public:
    typedef R ResultType;
    typedef typename TypeTraits<P1>::ParameterType Parm1;
    virtual R operator()(Parm1) = 0;
  };

////////////////////////////////////////////////////////////////////////////////
// class template FunctorImpl
// Specialization for 2 parameters
////////////////////////////////////////////////////////////////////////////////

  template <typename R, typename P1, typename P2>
  class FunctorImpl<R, TYPELIST_2(P1, P2)>
    : public Private::FunctorImplBase<R>
  {
  public:
    typedef R ResultType;
    typedef typename TypeTraits<P1>::ParameterType Parm1;
    typedef typename TypeTraits<P2>::ParameterType Parm2;
    virtual R operator()(Parm1, Parm2) = 0;
  };

////////////////////////////////////////////////////////////////////////////////
// class template FunctorImpl
// Specialization for 3 parameters
////////////////////////////////////////////////////////////////////////////////

  template <typename R, typename P1, typename P2, typename P3>
  class FunctorImpl<R, TYPELIST_3(P1, P2, P3)>
    : public Private::FunctorImplBase<R>
  {
  public:
    typedef R ResultType;
    typedef typename TypeTraits<P1>::ParameterType Parm1;
    typedef typename TypeTraits<P2>::ParameterType Parm2;
    typedef typename TypeTraits<P3>::ParameterType Parm3;
    virtual R operator()(Parm1, Parm2, Parm3) = 0;
  };

////////////////////////////////////////////////////////////////////////////////
// class template FunctorImpl
// Specialization for 4 parameters
////////////////////////////////////////////////////////////////////////////////

  template <typename R, typename P1, typename P2, typename P3, typename P4>
  class FunctorImpl<R, TYPELIST_4(P1, P2, P3, P4)>
    : public Private::FunctorImplBase<R>
  {
  public:
    typedef R ResultType;
    typedef typename TypeTraits<P1>::ParameterType Parm1;
    typedef typename TypeTraits<P2>::ParameterType Parm2;
    typedef typename TypeTraits<P3>::ParameterType Parm3;
    typedef typename TypeTraits<P4>::ParameterType Parm4;
    virtual R operator()(Parm1, Parm2, Parm3, Parm4) = 0;
  };

////////////////////////////////////////////////////////////////////////////////
// class template FunctorImpl
// Specialization for 5 parameters
////////////////////////////////////////////////////////////////////////////////

  template <typename R, typename P1, typename P2, typename P3, typename P4,
    typename P5>
  class FunctorImpl<R, TYPELIST_5(P1, P2, P3, P4, P5)>
    : public Private::FunctorImplBase<R>
  {
  public:
    typedef R ResultType;
    typedef typename TypeTraits<P1>::ParameterType Parm1;
    typedef typename TypeTraits<P2>::ParameterType Parm2;
    typedef typename TypeTraits<P3>::ParameterType Parm3;
    typedef typename TypeTraits<P4>::ParameterType Parm4;
    typedef typename TypeTraits<P5>::ParameterType Parm5;
    virtual R operator()(Parm1, Parm2, Parm3, Parm4, Parm5) = 0;
  };


////////////////////////////////////////////////////////////////////////////////
// class template FunctorHandler
// Wraps functors and pointers to functions
////////////////////////////////////////////////////////////////////////////////

  template <class ParentFunctor, typename Fun>
  class FunctorHandler
    : public ParentFunctor::Impl
  {
    typedef typename ParentFunctor::Impl Base;

  public:
    typedef typename Base::ResultType ResultType;
    typedef typename Base::Parm1 Parm1;
    typedef typename Base::Parm2 Parm2;
    typedef typename Base::Parm3 Parm3;
    typedef typename Base::Parm4 Parm4;
    typedef typename Base::Parm5 Parm5;
        
    FunctorHandler(const Fun& fun) : f_(fun) {}
        
    DEFINE_CLONE_FUNCTORIMPL(FunctorHandler)

      // operator() implementations for up to 5 arguments
                
      ResultType operator()()
      { return f_(); }

    ResultType operator()(Parm1 p1)
      { return f_(p1); }
        
    ResultType operator()(Parm1 p1, Parm2 p2)
      { return f_(p1, p2); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3)
      { return f_(p1, p2, p3); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3, Parm4 p4)
      { return f_(p1, p2, p3, p4); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3, Parm4 p4, Parm5 p5)
      { return f_(p1, p2, p3, p4, p5); }
        
  private:
    Fun f_;
  };
        
////////////////////////////////////////////////////////////////////////////////
// class template FunctorHandler
// Wraps pointers to member functions
////////////////////////////////////////////////////////////////////////////////

  template <class ParentFunctor, typename PointerToObj,
    typename PointerToMemFn>
  class MemFunHandler : public ParentFunctor::Impl
  {
    typedef typename ParentFunctor::Impl Base;

  public:
    typedef typename Base::ResultType ResultType;
    typedef typename Base::Parm1 Parm1;
    typedef typename Base::Parm2 Parm2;
    typedef typename Base::Parm3 Parm3;
    typedef typename Base::Parm4 Parm4;
    typedef typename Base::Parm5 Parm5;

    MemFunHandler(const PointerToObj& pObj, PointerToMemFn pMemFn) 
      : pObj_(pObj), pMemFn_(pMemFn)
      {}
        
    DEFINE_CLONE_FUNCTORIMPL(MemFunHandler)
        
      ResultType operator()()
      { return ((*pObj_).*pMemFn_)(); }

    ResultType operator()(Parm1 p1)
      { return ((*pObj_).*pMemFn_)(p1); }
        
    ResultType operator()(Parm1 p1, Parm2 p2)
      { return ((*pObj_).*pMemFn_)(p1, p2); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3)
      { return ((*pObj_).*pMemFn_)(p1, p2, p3); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3, Parm4 p4)
      { return ((*pObj_).*pMemFn_)(p1, p2, p3, p4); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3, Parm4 p4, Parm5 p5)
      { return ((*pObj_).*pMemFn_)(p1, p2, p3, p4, p5); }
        
  private:
    PointerToObj pObj_;
    PointerToMemFn pMemFn_;
  };
        
////////////////////////////////////////////////////////////////////////////////
// class template Functor
// A generalized functor implementation with value semantics
////////////////////////////////////////////////////////////////////////////////
       
  template <typename R, class TList = NullType>
  class ObjectFunctor
  {
  public:
    // Handy type definitions for the body type
    typedef FunctorImpl<R, TList> Impl;
    typedef R ResultType;
    typedef TList ParmList;
    typedef typename Impl::Parm1 Parm1;
    typedef typename Impl::Parm2 Parm2;
    typedef typename Impl::Parm3 Parm3;
    typedef typename Impl::Parm4 Parm4;
    typedef typename Impl::Parm5 Parm5;

    // Member functions

    ObjectFunctor() : spImpl_(0)
      {}
        
    ObjectFunctor(const ObjectFunctor& rhs) : spImpl_(Impl::Clone(rhs.spImpl_.get()))
      {}
        
    ObjectFunctor(std::auto_ptr<Impl> spImpl) : spImpl_(spImpl)
      {}
        
    template <typename Fun>
    ObjectFunctor(Fun fun)
      : spImpl_(new FunctorHandler<ObjectFunctor, Fun>(fun))
      {}

    template <class PtrObj, typename MemFn>
    ObjectFunctor(const PtrObj& p, MemFn memFn)
      : spImpl_(new MemFunHandler<ObjectFunctor, PtrObj, MemFn>(p, memFn))
      {}

    ObjectFunctor& operator=(const ObjectFunctor& rhs)
      {
	ObjectFunctor copy(rhs);
	// swap auto_ptrs by hand
	Impl* p = spImpl_.release();
	spImpl_.reset(copy.spImpl_.release());
	copy.spImpl_.reset(p);
	return *this;
      }
        
    ResultType operator()()
      { return (*spImpl_)(); }

    ResultType operator()(Parm1 p1)
      { return (*spImpl_)(p1); }
        
    ResultType operator()(Parm1 p1, Parm2 p2)
      { return (*spImpl_)(p1, p2); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3)
      { return (*spImpl_)(p1, p2, p3); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3, Parm4 p4)
      { return (*spImpl_)(p1, p2, p3, p4); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3, Parm4 p4, Parm5 p5)
      { return (*spImpl_)(p1, p2, p3, p4, p5); }

  private:
    std::auto_ptr<Impl> spImpl_;
  };
    
  namespace Private
  {
    template <class Fctor> struct BinderFirstTraits;

    template <typename R, class TList>
    struct BinderFirstTraits< ObjectFunctor<R, TList> >
    {
      typedef typename TL::Erase<TList, 
	typename TL::TypeAt<TList, 0>::Result>::Result
      ParmList;
      typedef ObjectFunctor<R, ParmList> BoundFunctorType;
      typedef typename BoundFunctorType::Impl Impl;
    };        
  }

////////////////////////////////////////////////////////////////////////////////
// class template BinderFirst
// Binds the first parameter of a ObjectFunctor object to a specific value
////////////////////////////////////////////////////////////////////////////////

  template <class OriginalFunctor>
  class BinderFirst 
    : public Private::BinderFirstTraits<OriginalFunctor>::Impl
  {
    typedef typename Private::BinderFirstTraits<OriginalFunctor>::Impl Base;
    typedef typename OriginalFunctor::ResultType ResultType;

    typedef typename OriginalFunctor::Parm1 BoundType;

    typedef typename OriginalFunctor::Parm2 Parm1;
    typedef typename OriginalFunctor::Parm3 Parm2;
    typedef typename OriginalFunctor::Parm4 Parm3;
    typedef typename OriginalFunctor::Parm5 Parm4;
    typedef EmptyType Parm5;

  public:
    BinderFirst(const OriginalFunctor& fun, BoundType bound)
      : f_(fun), b_(bound)
      {}
        
    DEFINE_CLONE_FUNCTORIMPL(BinderFirst)
        
      // operator() implementations for up to 5 arguments
                
      ResultType operator()()
      { return f_(b_); }

    ResultType operator()(Parm1 p1)
      { return f_(b_, p1); }
        
    ResultType operator()(Parm1 p1, Parm2 p2)
      { return f_(b_, p1, p2); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3)
      { return f_(b_, p1, p2, p3); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3, Parm4 p4)
      { return f_(b_, p1, p2, p3, p4); }
       
  private:
    OriginalFunctor f_;
    BoundType b_;
  };
    
////////////////////////////////////////////////////////////////////////////////
// function template BindFirst
// Binds the first parameter of a Functor object to a specific value
////////////////////////////////////////////////////////////////////////////////

  template <class Fctor>
  typename Private::BinderFirstTraits<Fctor>::BoundFunctorType
  BindFirst(
    const Fctor& fun,
    typename Fctor::Parm1 bound)
  {
    typedef typename Private::BinderFirstTraits<Fctor>::BoundFunctorType
      Outgoing;
        
    return Outgoing(std::auto_ptr<typename Outgoing::Impl>(
      new BinderFirst<Fctor>(fun, bound)));
  }

////////////////////////////////////////////////////////////////////////////////
// class template Chainer
// Chains two functor calls one after another
////////////////////////////////////////////////////////////////////////////////

  template <typename Fun1, typename Fun2>
  class Chainer : public Fun2::Impl
  {
    typedef Fun2 Base;

  public:
    typedef typename Base::ResultType ResultType;
    typedef typename Base::Parm1 Parm1;
    typedef typename Base::Parm2 Parm2;
    typedef typename Base::Parm3 Parm3;
    typedef typename Base::Parm4 Parm4;
    typedef typename Base::Parm5 Parm5;
        
    Chainer(const Fun1& fun1, const Fun2& fun2) : f1_(fun1), f2_(fun2) {}

    DEFINE_CLONE_FUNCTORIMPL(Chainer)

      // operator() implementations for up to 5 arguments

      ResultType operator()()
      { return f1_(), f2_(); }

    ResultType operator()(Parm1 p1)
      { return f1_(p1), f2_(p1); }
        
    ResultType operator()(Parm1 p1, Parm2 p2)
      { return f1_(p1, p2), f2_(p1, p2); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3)
      { return f1_(p1, p2, p3), f2_(p1, p2, p3); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3, Parm4 p4)
      { return f1_(p1, p2, p3, p4), f2_(p1, p2, p3, p4); }
        
    ResultType operator()(Parm1 p1, Parm2 p2, Parm3 p3, Parm4 p4, Parm5 p5)
      { return f1_(p1, p2, p3, p4, p5), f2_(p1, p2, p3, p4, p5); }
        
  private:
    Fun1 f1_;
    Fun2 f2_;
  };
    
////////////////////////////////////////////////////////////////////////////////
// function template Chain
// Chains two functor calls one after another
////////////////////////////////////////////////////////////////////////////////


  template <class Fun1, class Fun2>
  Fun2 Chain(
    const Fun1& fun1,
    const Fun2& fun2)
  {
    return Fun2(std::auto_ptr<typename Fun2::Impl>(
      new Chainer<Fun1, Fun2>(fun1, fun2)));
  }

}  // namespace Chroma

#endif
