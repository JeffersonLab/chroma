// -*- C++ -*-
// $Id: map_obj.h,v 1.1 2008-07-21 02:32:24 edwards Exp $
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
  private:
    //! Map type convenience
    typedef std::map<K,V> MapType_t;

  public:
    //! Default constructor
    MapObject() {}

    //! Destructor
    ~MapObject() {}

    //! Exists?
    bool exist(const K& key) const {
      return (src_map.find(key) == src_map.end()) ? false : true;
    }
			
    //! Insert
    void insert(const K& key, const V& val) {
      src_map.insert(std::make_pair(key,val));
    }
			
    //! Accessor
    const V& operator[](const K& key) const {
      if (! exist(key) )
      {
	QDPIO::cerr << "MapObject: key not found" << std::endl;
// No generic key writer
//	QDPIO::cerr << "key= " << key << std::endl;
//	QDPIO::cerr << "All Keys:" << std::endl;
//	std::vector<K> all_keys = dump();
//	for(int i=0; i < all_keys.size(); ++i)
//	  QDPIO::cerr << all_keys[i];

	exit(1);
      }

      return src_map.find(key)->second;
    }
			
    //! Dump keys
    std::vector<K> dump() const {
      std::vector<K> keys;
      typename MapType_t::const_iterator iter;
      for(iter  = src_map.begin();
	  iter != src_map.end();
	  ++iter)
      {
	keys.push_back(iter->first);
      }
      return keys;
    }

  private:
    //! Map of objects
    mutable MapType_t src_map;
  };

} // namespace Chroma

#endif
