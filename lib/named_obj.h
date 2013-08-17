// -*- C++ -*-
// $Id: named_obj.h,v 3.1 2007-09-22 01:19:03 edwards Exp $

/*! @file
 * @brief Named object support
 */

/*! \defgroup support Support routines
 * \ingroup lib
 *
 * Support routines
 */

#ifndef __named_obj_h__
#define __named_obj_h__

#include "chromabase.h"
#include "handle.h"
#include <map>
#include <string>

namespace Chroma
{
  //--------------------------------------------------------------------------------------
  //! Typeinfo Hiding Base Clase
  /*! @ingroup support
   */
  class NamedObjectBase 
  {
  public:
    NamedObjectBase() {}

    //! Setter
    virtual void setFileXML(XMLReader& xml) = 0;

    //! Setter
    virtual void setFileXML(XMLBufferWriter& xml) = 0;

    //! Setter
    virtual void setRecordXML(XMLReader& xml) = 0;

    //! Setter
    virtual void setRecordXML(XMLBufferWriter& xml) = 0;

    //! Getter
    virtual void getFileXML(XMLReader& xml) const = 0;

    //! Getter
    virtual void getFileXML(XMLBufferWriter& xml) const = 0;

    //! Getter
    virtual void getRecordXML(XMLReader& xml) const = 0;

    //! Getter
    virtual void getRecordXML(XMLBufferWriter& xml) const = 0;

    // This is key for cleanup
    virtual ~NamedObjectBase() {}
  };


  //--------------------------------------------------------------------------------------
  //! Type specific named object
  /*! @ingroup support
   */
  template<typename T>
  class NamedObject : public NamedObjectBase 
  {
  public:
    //! Constructor
    NamedObject() : data(new T) {}
  
    template<typename P1>
    NamedObject(const P1& p1) : data(new T(p1)) {}
 
    //! Destructor
    ~NamedObject() {}

    //! Setter
    void setFileXML(XMLReader& xml) 
    { 
      std::ostringstream os;
      xml.printCurrentContext(os);
      file_xml = os.str();
    }

    //! Setter
    void setFileXML(XMLBufferWriter& xml) 
    {
      file_xml = xml.printCurrentContext();
    }

    //! Setter
    void setRecordXML(XMLReader& xml) 
    {
      std::ostringstream os;
      xml.printCurrentContext(os);
      record_xml = os.str();
    }

    //! Setter
    void setRecordXML(XMLBufferWriter& xml) 
    {
      record_xml = xml.printCurrentContext();
    }

    //! Getter
    void getFileXML(XMLReader& xml) const 
    {
      std::istringstream os(file_xml);
      xml.open(os);
    }

    //! Getter
    void getFileXML(XMLBufferWriter& xml) const 
    {
      xml.writeXML(file_xml);
    }

    //! Getter
    void getRecordXML(XMLReader& xml) const 
    {
      std::istringstream os(record_xml);
      xml.open(os);
    }

    //! Getter
    void getRecordXML(XMLBufferWriter& xml) const 
    {
      xml.writeXML(record_xml);
    }

    //! Mutable data ref
    virtual T& getData() {
      return *data;
    }

    //! Const data ref
    virtual const T& getData() const {
      return *data;
    }

  private:
    Handle<T>   data;
    std::string file_xml;
    std::string record_xml;
  };


  //--------------------------------------------------------------------------------------
  //! The Map Itself
  /*! @ingroup support
   */
  class NamedObjectMap 
  {
  public:
    // Creation: clear the map
    NamedObjectMap() {
      the_map.clear();
    };

    // Destruction: erase all elements of the map
    ~NamedObjectMap() 
    {
      typedef std::map<std::string, NamedObjectBase*>::iterator I;
      while( ! the_map.empty() ) 
      {
	I iter = the_map.begin();

	delete iter->second;

	the_map.erase(iter);
      }
    }


    //! Create an entry of arbitrary type.
    template<typename T>
    void create(const std::string& id) 
    {
      // Lookup and throw exception if duplicate found
      typedef std::map<std::string, NamedObjectBase*>::iterator I;
      I iter = the_map.find(id);
      if( iter != the_map.end()) 
      {
	ostringstream error_stream;
        error_stream << "NamedObjectMap::create : duplicate id = " << id << endl;
        throw error_stream.str();
      }

      // Create a new object of specified type (empty)
      // Dynamic cast to Typeless base clasee.
      // Note multi1d's need to be loked up and resized appropriately
      // and no XML files are added at this point
      the_map[id] = dynamic_cast<NamedObjectBase*>(new NamedObject<T>());
      if (NULL == the_map[id])
      {
	ostringstream error_stream;
        error_stream << "NamedObjectMap::create : error creating NamedObject for id= " << id << endl;
        throw error_stream.str();
      }
    }

    //! Create an entry of arbitrary type, with 1 parameter
    template<typename T, typename P1>
    void create(const std::string& id, const P1& p1) 
    {
      // Lookup and throw exception if duplicate found
      MapType_t::iterator iter = the_map.find(id);
      if(iter != the_map.end()) 
      {
	ostringstream error_stream;
        error_stream << "NamedObjectMap::create : duplicate id = " << id << endl;
        throw error_stream.str();
      }

      // Create a new object of specified type (empty)
      // Dynamic cast to Typeless base clasee.
      // Note multi1d's need to be loked up and resized appropriately
      // and no XML files are added at this point
      the_map[id] = dynamic_cast<NamedObjectBase*>(new NamedObject<T>(p1));
      if (NULL == the_map[id])
      {
	ostringstream error_stream;
        error_stream << "NamedObjectMap::create : error creating NamedObject for id= " << id << endl;
        throw error_stream.str();
      }
    }


    //! Check if an id exists
    bool check(const std::string& id) const
    {
      // Do a lookup
      MapType_t::const_iterator iter = the_map.find(id);
    
      // If found then return true
      return (iter != the_map.end()) ? true : false;
    }
  
  
    //! Delete an item that we no longer neeed
    void erase(const std::string& id) 
    {
      // Do a lookup
      MapType_t::iterator iter = the_map.find(id);
    
      // If found then delete it.
      if( iter != the_map.end() ) 
      { 
      	// Delete the data.of the record
	delete iter->second;

	// Delete the record
	the_map.erase(iter);
      }
      else 
      {
	// We attempt to erase something non existent
	ostringstream error_stream;
        error_stream << "NamedObjectMap::erase : erasing unknown id = " << id << endl;
        throw error_stream.str();
      }
    }
  
  
    //! Dump out all objects
    void dump() const
    {
      QDPIO::cout << "Available Keys are : " << endl;
      for(MapType_t::const_iterator j = the_map.begin(); j != the_map.end(); j++) 
	QDPIO::cout << j->first << endl;
    }
  
  
    //! Look something up and return a NamedObjectBase reference
    NamedObjectBase& get(const std::string& id) const
    {
      // Find it
      MapType_t::const_iterator iter = the_map.find(id);
      if (iter == the_map.end()) 
      {
	// Not found -- lookup exception
	ostringstream error_stream;
        error_stream << "NamedObjectMap::get : unknown id = " << id << endl;
        throw error_stream.str();
      }
      else 
      {
	// Found, return the reference
	return *(iter->second);
      }
    }

    //! Look something up and return a ref to the derived named object
    template<typename T>
    T& getObj(const std::string& id) 
    {
      return dynamic_cast<NamedObject<T>&>(get(id));
    }

    //! Look something up and return a ref to the derived named object
    template<typename T>
    const T& getObj(const std::string& id) const 
    {
      return dynamic_cast<NamedObject<T>&>(get(id));
    }

    //! Look something up and return a ref to actual data
    template<typename T>
    T& getData(const std::string& id) 
    {
      return dynamic_cast<NamedObject<T>&>(get(id)).getData();
    }

    //! Look something up and return a ref to actual data
    template<typename T>
    const T& getData(const std::string& id) const 
    {
      return dynamic_cast<NamedObject<T>&>(get(id)).getData();
    }

  private:
    typedef std::map<std::string, NamedObjectBase*> MapType_t;
    MapType_t the_map;
  };

}

#endif
