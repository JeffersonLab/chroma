// -*- C++ -*-
/*! \file map_obj_memory.h
 *  \brief A memory based map object 
 */

#ifndef __map_obj_memory_h__
#define __map_obj_memory_h__

#include "util/ferm/map_obj.h"

namespace Chroma
{

  //! Private Namespace 
  namespace MapObjectMemoryEnv { 
    
    //! Register the MapObject in factories
    bool registerAll();

  }

  //----------------------------------------------------------------------------
  //! A wrapper over maps
  template<typename K, typename V>
  class MapObjectMemory : public MapObject<K,V>
  {
  public:
    //! Map type convenience
    typedef std::map<K,V> MapType_t;

    //! Default constructor
    MapObjectMemory() {}

    //! Destructor
    ~MapObjectMemory() {}

    //! Exists?
    bool exist(const K& key) const {
      return (src_map.find(key) == src_map.end()) ? false : true;
    }
			
    //! Insert
    void insert(const K& key, const V& val) {
      src_map.insert(std::make_pair(key,val));
    }
			
    //! Accessor
    void lookup(const K& key, V& val) const { 
      if (! exist(key) ) {
	QDPIO::cerr << "MapObjectMemory: key not found" << std::endl;
// No generic key writer
//	QDPIO::cerr << "key= " << key << std::endl;
//	QDPIO::cerr << "All Keys:" << std::endl;
//	std::vector<K> all_keys = dump();
//	for(int i=0; i < all_keys.size(); ++i)
//	  QDPIO::cerr << all_keys[i];

	exit(1);
      }

      val = src_map.find(key)->second;
    }
			
    //! The number of elements
    typename MapType_t::size_type size() const {return src_map.size();}

    //! Dump keys
    std::vector<K> dump() const {
      std::vector<K> keys;
      typename MapType_t::const_iterator iter;
      for(iter  = src_map.begin();
	  iter != src_map.end();
	  ++iter) { 
	keys.push_back(iter->first);
      }
      return keys;
    }
    
    //! Usual begin iterator
    typename MapType_t::const_iterator begin() const {return src_map.begin();}

    //! Usual end iterator
    typename MapType_t::const_iterator end() const {return src_map.end();}

  private:
    //! Map of objects
    mutable MapType_t src_map;
  };

} // namespace Chroma

#endif
