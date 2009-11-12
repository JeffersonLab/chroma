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

    //! Default constructor
    MapObject() {}

    //! Virtual Destructor
    virtual	
    ~MapObject() {}

    //! Open object in read mode (lookups)
    virtual
    bool openRead(void) = 0;
   
    //! Close read mode
    virtual 
    bool closeRead(void) = 0;

    
    //! Open write mode (inserts)
    virtual
    bool openWrite(void) = 0;
    
    //! Close write mode
    virtual
    bool closeWrite(void) = 0;


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
    unsigned long size() const = 0;

    //! Dump keys
    virtual
    std::vector<K> dump() const = 0;

  };

} // namespace Chroma

#endif
