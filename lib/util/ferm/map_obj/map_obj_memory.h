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
    MapObjectMemory() : readMode(false), writeMode(false), state(INIT) {}


    //! open Read mode (Lookups)
    void openRead(void) {  
      state = READ;
    }


    //! open Update mode (updates)
    void openUpdate(void) {
      switch(state) { 
      case UPDATE:// Deliberate fall through
      case READ:
	state = UPDATE;
	break;
      case INIT:
	errorState("MapObjectMemory: openUpdate() cannot update from INIT mode. Write first");
	break;
      case WRITE:
	errorState("MapObjectMemory: openUpdate() cannot update from WRITE mode. Finish writing first");
	break;
      default:
	errorState("MapObjectMemory: openUpdate() called from unknown mode");
	break;
      }
    }
      

    //! open Write mode (Inserts)
    void openWrite(void) {
      switch(state) {
      case INIT:
	state = WRITE;
	break;
      default:
	errorState("MapObjectMemory: openWrite() should only be called "
		   "from INIT mode");
	break;
      }
    }

    //! Destructor
    ~MapObjectMemory() {}


    //! Insert
    void insert(const K& key, const V& val) {
      if(state==WRITE) {
	src_map.insert(std::make_pair(key,val));
      }
      else { 
	errorState("MapObjectMemory: insert() called from outsite write mode");
      }
    }

    //! Update
    void update(const K& key, const V& val) { 
      if(state==UPDATE) {
	if( exist(key) ) {
	  src_map[key] = val; // Update
	}
	else { 
	  errorState("MapObjectMemory: your key is not in the map. You cannot update a non existent key");
	}
      }
      else { 
	errorState("MapObjectMemory: update() can only be called from Update mode");
      }
    }

    //! Accessor
    void lookup(const K& key, V& val) const { 
      if(state == READ || state==UPDATE) { 
	if (! exist(key) ) {
	  QDPIO::cout << "Couldnt find key " <<endl;
	  dump();
	  errorState("MapObjectMemory: lookup: key_not_found"); 
	}
	
	val = src_map.find(key)->second;
      }
      else { 
	errorState("MapObjectMemory: lookup called from outside READ mode");
      }
    }

    
    //! Exists?
    bool exist(const K& key) const {
      return (src_map.find(key) == src_map.end()) ? false : true;
    }
			
    //! The number of elements
    unsigned int size() const {return static_cast<unsigned long>(src_map.size());}

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

    //! An error sink state. Throws exception
    void errorState(const std::string& e) const {
      throw e;
    }

  private:
    //! Map of objects
    mutable MapType_t src_map;

    mutable bool readMode;
    mutable bool writeMode;

    enum State { INIT, READ, WRITE, UPDATE};
    State state;

  };

} // namespace Chroma

#endif
