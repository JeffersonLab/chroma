// -*- C++ -*-
/*! @file
 * @brief Singleton support
 */

#ifndef __singleton_h__
#define __singleton_h__

#include "qdp.h"

#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <cstdlib>
#include <new>

namespace Chroma
{

////////////////////////////////////////////////////////////////////////////////
// class template SingleThreaded
// Implementation of the ThreadingModel policy used by various classes
// Implements a single-threaded model; no synchronization
////////////////////////////////////////////////////////////////////////////////

  template <class Host>
  class SingleThreaded
  {
  public:
    struct Lock
    {
      Lock() {}
      Lock(const Host&) {}
    };
        
    typedef Host VolatileType;

    typedef int IntType; 

    static IntType AtomicAdd(volatile IntType& lval, IntType val)
      { return lval += val; }
        
    static IntType AtomicSubtract(volatile IntType& lval, IntType val)
      { return lval -= val; }

    static IntType AtomicMultiply(volatile IntType& lval, IntType val)
      { return lval *= val; }
        
    static IntType AtomicDivide(volatile IntType& lval, IntType val)
      { return lval /= val; }
        
    static IntType AtomicIncrement(volatile IntType& lval)
      { return ++lval; }
        
    static IntType AtomicDivide(volatile IntType& lval)
      { return --lval; }
        
    static void AtomicAssign(volatile IntType & lval, IntType val)
      { lval = val; }
        
    static void AtomicAssign(IntType & lval, volatile IntType & val)
      { lval = val; }
  };



  namespace Private
  {
////////////////////////////////////////////////////////////////////////////////
// class LifetimeTracker
// Helper class for SetLongevity
////////////////////////////////////////////////////////////////////////////////

    class LifetimeTracker
    {
    public:
      LifetimeTracker(unsigned int x) : longevity_(x) 
	{}
            
      virtual ~LifetimeTracker() = 0;
            
      static bool Compare(const LifetimeTracker* lhs,
			  const LifetimeTracker* rhs)
	{
	  return rhs->longevity_ > lhs->longevity_;
	}
            
    private:
      unsigned int longevity_;
    };
        
    // Definition required
    inline LifetimeTracker::~LifetimeTracker() {} 
        
    // Helper data
    typedef LifetimeTracker** TrackerArray;
    extern TrackerArray pTrackerArray;
    extern unsigned int elements;

    // Helper destroyer function
    template <typename T>
    struct Deleter
    {
      static void Delete(T* pObj)
	{ delete pObj; }
    };

    // Concrete lifetime tracker for objects of type T
    template <typename T, typename Destroyer>
    class ConcreteLifetimeTracker : public LifetimeTracker
    {
    public:
      ConcreteLifetimeTracker(T* p,unsigned int longevity, Destroyer d)
	: LifetimeTracker(longevity)
	, pTracked_(p)
	, destroyer_(d)
	{}
            
      ~ConcreteLifetimeTracker()
	{ destroyer_(pTracked_); }
            
    private:
      T* pTracked_;
      Destroyer destroyer_;
    };

    void AtExitFn(); // declaration needed below
    
  } // namespace Private

////////////////////////////////////////////////////////////////////////////////
// function template SetLongevity
// Assigns an object a longevity; ensures ordered destructions of objects 
//     registered thusly during the exit sequence of the application
////////////////////////////////////////////////////////////////////////////////

  template <typename T, typename Destroyer>
  void SetLongevity(T* pDynObject, unsigned int longevity,
		    Destroyer d = Private::Deleter<T>::Delete)
  {
    using namespace Private;
        
    TrackerArray pNewArray = static_cast<TrackerArray>(
      std::realloc(pTrackerArray, elements + 1));
    if (!pNewArray) throw std::bad_alloc();
        
    LifetimeTracker* p = new ConcreteLifetimeTracker<T, Destroyer>(
      pDynObject, longevity, d);
        
    // Delayed assignment for exception safety
    pTrackerArray = pNewArray;
        
    // Insert a pointer to the object into the queue
    TrackerArray pos = std::upper_bound(
      pTrackerArray, 
      pTrackerArray + elements, 
      p, 
      LifetimeTracker::Compare);
    std::copy_backward(
      pos, 
      pTrackerArray + elements,
      pTrackerArray + elements + 1);
    *pos = p;
    ++elements;
        
    // Register a call to AtExitFn
    std::atexit(Private::AtExitFn);
  }

////////////////////////////////////////////////////////////////////////////////
// class template CreateUsingNew
// Implementation of the CreationPolicy used by SingletonHolder
// Creates objects using a straight call to the new operator 
////////////////////////////////////////////////////////////////////////////////

  template <class T> struct CreateUsingNew
  {
    static T* Create()
      { return new T; }
        
    static void Destroy(T* p)
      { delete p; }
  };
    
////////////////////////////////////////////////////////////////////////////////
// class template CreateUsingNew
// Implementation of the CreationPolicy used by SingletonHolder
// Creates objects using a call to std::malloc, followed by a call to the 
//     placement new operator
////////////////////////////////////////////////////////////////////////////////

  template <class T> struct CreateUsingMalloc
  {
    static T* Create()
      {
	void* p = std::malloc(sizeof(T));
	if (!p) return 0;
	return new(p) T;
      }
        
    static void Destroy(T* p)
      {
	p->~T();
	std::free(p);
      }
  };
    
////////////////////////////////////////////////////////////////////////////////
// class template CreateStatic
// Implementation of the CreationPolicy used by SingletonHolder
// Creates an object in static memory
// Implementation is slightly nonportable because it uses the MaxAlign trick 
//     (an union of all types to ensure proper memory alignment). This trick is 
//     nonportable in theory but highly portable in practice.
////////////////////////////////////////////////////////////////////////////////

  template <class T> struct CreateStatic
  {
    union MaxAlign
    {
      char t_[sizeof(T)];
      short int shortInt_;
      int int_;
      long int longInt_;
      float float_;
      double double_;
      long double longDouble_;
      struct Test;
      int Test::* pMember_;
      int (Test::*pMemberFn_)(int);
    };
        
    static T* Create()
      {
	static MaxAlign staticMemory_;
	return new(&staticMemory_) T;
      }
        
    static void Destroy(T* p)
      {
	p->~T();
      }
  };
    
////////////////////////////////////////////////////////////////////////////////
// class template DefaultLifetime
// Implementation of the LifetimePolicy used by SingletonHolder
// Schedules an object's destruction as per C++ rules
// Forwards to std::atexit
////////////////////////////////////////////////////////////////////////////////

  template <class T>
  struct DefaultLifetime
  {
    static void ScheduleDestruction(T*, void (*pFun)())
      { std::atexit(pFun); }
        
    static void OnDeadReference()
      { throw std::logic_error("Dead Reference Detected"); }
  };

  // Copy to help with disambiguation
  template <class T>
  struct DefaultLifetime1
  {
    static void ScheduleDestruction(T*, void (*pFun)())
      { std::atexit(pFun); }
        
    static void OnDeadReference()
      { throw std::logic_error("Dead Reference Detected"); }
  };

  // Copy to help with disambiguation
  template <class T>
  struct DefaultLifetime2
  {
    static void ScheduleDestruction(T*, void (*pFun)())
      { std::atexit(pFun); }
        
    static void OnDeadReference()
      { throw std::logic_error("Dead Reference Detected"); }
  };

////////////////////////////////////////////////////////////////////////////////
// class template PhoenixSingleton
// Implementation of the LifetimePolicy used by SingletonHolder
// Schedules an object's destruction as per C++ rules, and it allows object 
//    recreation by not throwing an exception from OnDeadReference
////////////////////////////////////////////////////////////////////////////////

  template <class T>
  class PhoenixSingleton
  {
  public:
    static void ScheduleDestruction(T*, void (*pFun)())
      {
#ifndef ATEXIT_FIXED
	if (!destroyedOnce_)
#endif
	  std::atexit(pFun);
      }
        
    static void OnDeadReference()
      {
#ifndef ATEXIT_FIXED
	destroyedOnce_ = true;
#endif
      }
        
  private:
#ifndef ATEXIT_FIXED
    static bool destroyedOnce_;
#endif
  };
    
#ifndef ATEXIT_FIXED
  template <class T> bool PhoenixSingleton<T>::destroyedOnce_ = false;
#endif
        
////////////////////////////////////////////////////////////////////////////////
// class template Adapter
// Helper for SingletonWithLongevity below
////////////////////////////////////////////////////////////////////////////////

  namespace Private
  {
    template <class T>
    struct Adapter
    {
      void operator()(T*) { return pFun_(); }
      void (*pFun_)();
    };
  }

////////////////////////////////////////////////////////////////////////////////
// class template SingletonWithLongevity
// Implementation of the LifetimePolicy used by SingletonHolder
// Schedules an object's destruction in order of their longevities
// Assumes a visible function GetLongevity(T*) that returns the longevity of the
//     object
////////////////////////////////////////////////////////////////////////////////

  template <class T>
  class SingletonWithLongevity
  {
  public:
    static void ScheduleDestruction(T* pObj, void (*pFun)())
      {
	Private::Adapter<T> adapter = { pFun };
	SetLongevity(pObj, GetLongevity(pObj), adapter);
      }
        
    static void OnDeadReference()
      { throw std::logic_error("Dead Reference Detected"); }
  };

////////////////////////////////////////////////////////////////////////////////
// class template NoDestroy
// Implementation of the LifetimePolicy used by SingletonHolder
// Never destroys the object
////////////////////////////////////////////////////////////////////////////////

  template <class T>
  struct NoDestroy
  {
    static void ScheduleDestruction(T*, void (*)())
      {}
        
    static void OnDeadReference()
      {}
  };

////////////////////////////////////////////////////////////////////////////////
// class template SingletonHolder
// Provides Singleton amenities for a type T
// To protect that type from spurious instantiations, you have to protect it
//     yourself.
////////////////////////////////////////////////////////////////////////////////

  template <typename T,
            template <class> class CreationPolicy = CreateUsingNew,
            template <class> class LifetimePolicy = DefaultLifetime,
            template <class> class ThreadingModel = SingleThreaded>
  class SingletonHolder
  {
  public:
    static T& Instance();
        
  private:
    // Helpers
    static void MakeInstance();
    static void DestroySingleton();
        
    // Protection
    SingletonHolder();
        
    // Data
    typedef typename ThreadingModel<T*>::VolatileType PtrInstanceType;
    static PtrInstanceType pInstance_;
    static bool destroyed_;
  };
    
////////////////////////////////////////////////////////////////////////////////
// SingletonHolder's data
////////////////////////////////////////////////////////////////////////////////

  template <class T,
            template <class> class C,
            template <class> class L,
            template <class> class M>
  typename SingletonHolder<T, C, L, M>::PtrInstanceType
  SingletonHolder<T, C, L, M>::pInstance_;

  template
  <
  class T,
    template <class> class C,
    template <class> class L,
    template <class> class M
  >
  bool SingletonHolder<T, C, L, M>::destroyed_;

////////////////////////////////////////////////////////////////////////////////
// SingletonHolder::Instance
////////////////////////////////////////////////////////////////////////////////

  template <class T,
            template <class> class CreationPolicy,
            template <class> class LifetimePolicy,
            template <class> class ThreadingModel>
  inline T& SingletonHolder<T, CreationPolicy, 
    LifetimePolicy, ThreadingModel>::Instance()
  {
    if (!pInstance_)
    {
      MakeInstance();
    }
    return *pInstance_;
  }

////////////////////////////////////////////////////////////////////////////////
// SingletonHolder::MakeInstance (helper for Instance)
////////////////////////////////////////////////////////////////////////////////

  template <class T,
            template <class> class CreationPolicy,
            template <class> class LifetimePolicy,
            template <class> class ThreadingModel>
  void SingletonHolder<T, CreationPolicy, 
    LifetimePolicy, ThreadingModel>::MakeInstance()
  {
    typename ThreadingModel<T>::Lock guard;
    (void)guard;
        
    if (!pInstance_)
    {
      if (destroyed_)
      {
	LifetimePolicy<T>::OnDeadReference();
	destroyed_ = false;
      }
      pInstance_ = CreationPolicy<T>::Create();
      LifetimePolicy<T>::ScheduleDestruction(pInstance_, 
					     &DestroySingleton);
    }
  }

  template <class T,
            template <class> class CreationPolicy,
            template <class> class L,
            template <class> class M>
  void SingletonHolder<T, CreationPolicy, L, M>::DestroySingleton()
  {
    assert(!destroyed_);
    CreationPolicy<T>::Destroy(pInstance_);
    pInstance_ = 0;
    destroyed_ = true;
  }
} // namespace Chroma


#endif
