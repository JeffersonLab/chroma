// -*- C++ -*-
// $Id: disk_obj.h,v 1.1 2009-10-30 17:12:55 edwards Exp $
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
  class DiskObject
  {
  private:
    //! Disk type convenience
    typedef std::map<K,int> DiskType_t;

  public:
    //! Default constructor
    DiskObject(const std::string& filename) {}

    //! Destructor
    ~DiskObject() {}

    //! Exists?
    bool exist(const K& key) const {
      return (src_map.find(key) == src_map.end()) ? false : true;
    }
			
    //! Insert
    void insert(const K& key, const V& val) {
      src_map.insert(std::make_pair(key,val));
    }
			
    //! Accessor
    void get(const K& key, V& val) const {
      if (! exist(key) )
      {
	QDPIO::cerr << "DiskObject: key not found" << std::endl;
// No generic key writer
//	QDPIO::cerr << "key= " << key << std::endl;
//	QDPIO::cerr << "All Keys:" << std::endl;
//	std::vector<K> all_keys = dump();
//	for(int i=0; i < all_keys.size(); ++i)
//	  QDPIO::cerr << all_keys[i];

	exit(1);
      }

      int recno = src_map.find(key)->second;
      
    }
			
    //! Accessor
    V operator[](const K& key) const {
      V val;
      get(key, val);
      return val;
    }

    //! The number of elements
    typename DiskType_t::size_type size() const {return src_map.size();}

    //! Dump keys
    std::vector<K> keys() const {
      std::vector<K> keys;
      typename DiskType_t::const_iterator iter;
      for(iter  = src_map.begin();
	  iter != src_map.end();
	  ++iter)
      {
	keys.push_back(iter->first);
      }
      return keys;
    }

  private:
    //! Disk of objects
    mutable DiskType_t src_map;
    QDPFileReader qdp_in;
    QDPFileWriter qdp_out;
  };

} // namespace Chroma

#endif
