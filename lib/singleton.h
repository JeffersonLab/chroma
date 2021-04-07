// -*- C++ -*-
/*! @file
 * @brief Singleton support
 */

#ifndef __singleton_h__
#define __singleton_h__

#include <functional>
#include <vector>

namespace Chroma {

// List of singleton to destroy at chroma finalize
using DestroyList = std::vector<std::function<void()>>;
inline DestroyList &getDestroyList() {
  static DestroyList list;
  return list;
}

// Destroy all singleton created with SingletonHolder
inline void destroySingletons() {
  // Destroy in reverse order to the allocated one
  while (getDestroyList().size() > 0) {
    getDestroyList().back()();
    getDestroyList().pop_back();
  }
}

////////////////////////////////////////////////////////////////////////////////
// class template SingletonHolder
// Provides Singleton amenities for a type T; instances are destroy at
// function finalize. Create several singletons of the same type T by setting
// different types for InstanceT.
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename InstanceT=void> class SingletonHolder {
public:
  static T &Instance() {
    static T *instance = nullptr;
    if (!instance) {
      instance = new T;
      getDestroyList().push_back([=]() { delete instance; });
    }
    return *instance;
  }

private:
  // Protection
  SingletonHolder() {}
};

} // namespace Chroma

#endif
