// -*- C++ -*-
/*! \file map_obj_memory.h
 *  \brief A memory based map object 
 */

#ifndef __map_obj_memory_h__
#define __map_obj_memory_h__

#include "util/ferm/map_obj.h"

namespace Chroma
{


  //----------------------------------------------------------------------------
  //! A wrapper over maps
  template<typename K, typename V>
  class MapObjectMemory : public MapObject<K,V>
  {
  public:

    //! Default constructor
    MapObjectMemory() : readMode(false), writeMode(false) {}


    //! OpenRead mode (Inserts) 
    bool openRead(void) {  
      bool ret_val=true;

      if ( !writeMode ) {
	if( !readMode )  { 
	  readMode = true;
	  ret_val = true;
	}
	else { 
	  QDPIO::cerr << "MapObjectMemory: Already Open in read mode" << endl;
	  ret_val = true; // Re-opening in same mode is not an error
	}
	  
      }
      else { 
	QDPIO::cerr << "MapObjectMemory: Already  open in write mode" << endl;
	readMode = false;
	ret_val = false; // ReadOpening without close from write mode is an error
      }
      
      return ret_val;
    }

    bool closeRead(void) { 
      bool ret_val = true;
      if( readMode ) {
	readMode = false; 
	ret_val = true;
      }
      else {
	QDPIO::cerr << "MapObjectMemory: can\'t closeRead if not in read mode" << endl;
	ret_val = false; // 
      }

      return ret_val;
    }

    bool openWrite(void) {
      bool ret_val = true;
      if( !readMode  ) {
	if( !writeMode ) { 
	  writeMode = true;
	  ret_val = true;
	}
	else { 
	  QDPIO::cerr << "MapObjectMemory: already open in write mode" << endl;
	  ret_val = true; // ReOpening from same mode is not really an error
	}
      }
      else { 
	QDPIO::cerr << "MapObjectMemory: can\'t OpenWrite if already open in read mode" << endl;
	ret_val = false;
      }
      return ret_val;
    }

    bool closeWrite(void) { 
      bool ret_val = true;
      if( writeMode ) {
	writeMode = false;
	ret_val = true;
      }
      else { 
	QDPIO::cerr << "Can\'t close write if not open in write mode" << endl;
	ret_val = false;
      }
      return ret_val;
    }

       

    //! Destructor
    ~MapObjectMemory() {
      if( readMode) closeRead();
      if( writeMode) closeWrite();
    }

    //! Exists?
    bool exist(const K& key) const {
      return (src_map.find(key) == src_map.end()) ? false : true;
    }
			
    //! Insert
    void insert(const K& key, const V& val) {
      if(writeMode) {
	src_map.insert(std::make_pair(key,val));
      }
      else { 
	QDPIO::cerr << "MapObjectMemory: Cant insert if not open in writeMode" << endl;
	throw std::string("MapObjectMemory: Cant insert without being in write mode");
      }
    }
			
    //! Accessor
    void lookup(const K& key, V& val) const { 
      if( readMode ) {
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
      else { 
	QDPIO::cerr << "MapObjectMemory: Cant lookup if not open in readMode" << endl;
	throw std::string("MapObjectMemory: Cant lookup if not open in readMode");
      }
    }
			
    //! The number of elements
    unsigned long size() const {return static_cast<unsigned long>(src_map.size());}

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

    /*! 
     * These extend the bacis MapObject Interface. 
     * The iterators are used to QIO the object
     * Need to be public for now 
     */

    //! Usual begin iterator
    //! Map type convenience
    typedef std::map<K,V> MapType_t;
    
 
    typename MapType_t::const_iterator begin() const {return src_map.begin();}
    
    //! Usual end iterator
    typename MapType_t::const_iterator  end() const {return src_map.end();}


  private:
    //! Map of objects
    mutable MapType_t src_map;

    mutable bool readMode;
    mutable bool writeMode;



  };

} // namespace Chroma

#endif
