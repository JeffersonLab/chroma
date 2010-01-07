// -*- C++ -*-
/*! \file
 *  \brief A null map object 
 */

#ifndef __map_obj_null_h__
#define __map_obj_null_h__

#include "util/ferm/map_obj.h"

namespace Chroma
{


  //----------------------------------------------------------------------------
  //! A wrapper over maps
  template<typename K, typename V>
  class MapObjectNull : public MapObject<K,V>
  {
  public:

    //! Default constructor
    MapObjectNull() : readMode(false), writeMode(false), state(INIT) {}


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
	errorState("MapObjectNull: openUpdate() cannot update from INIT mode. Write first");
	break;
      case WRITE:
	errorState("MapObjectNull: openUpdate() cannot update from WRITE mode. Finish writing first");
	break;
      default:
	errorState("MapObjectNull: openUpdate() called from unknown mode");
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
	errorState("MapObjectNull: openWrite() should only be called from INIT mode");
	break;
      }
    }

    //! Destructor
    ~MapObjectNull() {}


    //! Insert
    void insert(const K& key, const V& val) {
      if(state==WRITE) {
	// throw away
      }
      else { 
	errorState("MapObjectNull: insert() called from outsite write mode");
      }
    }

    //! Update
    void update(const K& key, const V& val) { 
      if(state==UPDATE) {
	errorState("MapObjectNull: does not support update");
      }
      else { 
	errorState("MapObjectNull: update() can only be called from Update mode");
      }
    }

    //! Accessor
    void lookup(const K& key, V& val) const { 
      if(state == READ || state==UPDATE) { 
	errorState("MapObjectNull: does not support lookup"); 
      }
      else { 
	errorState("MapObjectNull: lookup called from outside READ mode");
      }
    }

    
    //! Exists?
    bool exist(const K& key) const {return false;}
			
    //! The number of elements
    unsigned int size() const {return 0;}

    //! Dump keys
    std::vector<K> dump() const {std::vector<K> keys; return keys;}

    /*! 
     * These extend the bacis MapObject Interface. 
     * The iterators are used to QIO the object
     * Need to be public for now 
     */

    //! Usual begin iterator
    //! Map type convenience
    typedef std::map<K,V> MapType_t;
    

    //! Annoying, need these to satisfy the map
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
