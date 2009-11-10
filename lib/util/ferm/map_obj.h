// -*- C++ -*-
/*! \file
 * \brief Wrapper over maps
 */

#ifndef __map_obj_h__
#define __map_obj_h__

#include <map>
#include <vector>

namespace Chroma
{

  //----------------------------------------------------------------------------
  //! A wrapper over maps
  template<typename K, typename V>
  class MapObject
  {
  public:
    //! Map type convenience
    typedef std::map<K,V> MapType_t;

    //! Default constructor
    MapObject() {}

    //! Virtual Destructor
    virtual	
    ~MapObject() {}

    //! Exists?
    virtual
    bool exist(const K& key) const = 0;
			
    //! Insert
    virtual
    void insert(const K& key, const V& val) = 0;
			
    //! Other accessor
    virtual
    void lookup(const K& key, V& val) const = 0;

    //! Size of Map
    virtual
    typename MapType_t::size_type size() const = 0;

    //! Dump keys
    virtual
    std::vector<K> dump() const = 0;

    //! Usual begin iterator
    virtual
    typename MapType_t::const_iterator begin() const = 0; 

    //! Usual end iterator
    virtual
    typename MapType_t::const_iterator end() const = 0;
  };

} // namespace Chroma

#endif
