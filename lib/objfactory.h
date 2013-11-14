// -*- C++ -*-
/*! @file
 * @brief Factory class for objects from XML input
 */

#ifndef __objfact_h__
#define __objfact_h__

#include "chromabase.h"

#include "typeinfo.h"
#include "typelist.h"
#include <map>
#include <exception>


namespace Chroma
{

////////////////////////////////////////////////////////////////////////////////
// class template DefaultFactoryError
// Manages the "Unknown Type" error in an object factory
////////////////////////////////////////////////////////////////////////////////

  template <typename IdentifierType, class AbstractProduct>
  struct DefaultFactoryError
  {
    struct Exception : public std::exception
    {
      const char* what() const throw() { return "Unknown Type"; }
    };
        
    static AbstractProduct* OnUnknownType(const IdentifierType& id)
      {
	throw Exception();
      }
  };

////////////////////////////////////////////////////////////////////////////////
// class template StringFactoryError
// Manages the "Unknown Type" error in an object factory
////////////////////////////////////////////////////////////////////////////////

  template <typename IdentifierType, class AbstractProduct>
  struct StringFactoryError
  {
    static AbstractProduct* OnUnknownType(const IdentifierType& id)
      {
	std::ostringstream os;
	os << "Factory error: unknown identifier: id = " << id << endl;
	throw os.str();
      }
  };

////////////////////////////////////////////////////////////////////////////////
// class template NullFactoryError
// Manages the "Unknown Type" error in an object factory
////////////////////////////////////////////////////////////////////////////////

  template <typename IdentifierType, class AbstractProduct>
  struct NullFactoryError
  {
    static AbstractProduct* OnUnknownType(const IdentifierType&)
      {
	return 0;
      }
  };


//! Object factory class
/*! @ingroup actions
 *
 * Supports abstract creation of objects
 */
  template<class AbstractProduct, 
           typename IdentifierType,
           class TList = NullType,
           typename ProductCreator = AbstractProduct* (*)(),
           template<typename, class>
             class FactoryErrorPolicy = DefaultFactoryError>
  class ObjectFactory : public FactoryErrorPolicy<IdentifierType, AbstractProduct>
  {
  public:
    // Handy type definitions for the body type
    typedef TList ParmList;
    typedef typename TL::TypeAtNonStrict<ParmList,0,NullType>::Result  Parm1;
    typedef typename TL::TypeAtNonStrict<ParmList,1,NullType>::Result  Parm2;
    typedef typename TL::TypeAtNonStrict<ParmList,2,NullType>::Result  Parm3;
    typedef typename TL::TypeAtNonStrict<ParmList,3,NullType>::Result  Parm4;
    typedef typename TL::TypeAtNonStrict<ParmList,4,NullType>::Result  Parm5;
    typedef typename TL::TypeAtNonStrict<ParmList,5,NullType>::Result  Parm6;
    typedef typename TL::TypeAtNonStrict<ParmList,6,NullType>::Result  Parm7;
    typedef typename TL::TypeAtNonStrict<ParmList,7,NullType>::Result  Parm8;
    typedef typename TL::TypeAtNonStrict<ParmList,8,NullType>::Result  Parm9;

    // Member functions
    //! Register the object
    /*! 
     * \param id       object id
     * \param creator  the callback function to create the object
     * \return returns true if registration successful
     */
    bool registerObject(const IdentifierType& id, ProductCreator creator)
      {
	return associations_.insert(
	  typename IdToProductMap::value_type(id, creator)).second;
      }
  
    //! Unregister the object
    /*! 
     * \param id       object id
     * \return returns true if object name was registered before
     */
    bool unregisterObject(const IdentifierType& id)
      {
	return associations_.erase(id) == 1;
      }

    //! Create the object
    /*! 
     * \param id       object id
     * \return returns pointer to the object
     */
    AbstractProduct* createObject(const IdentifierType& id)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) {
	
	typedef typename IdToProductMap::const_iterator CI;
	QDPIO::cerr << "Couldnt find key " << id << " in the map: " << endl;
	QDPIO::cerr << "Available Keys are : " << endl;
	for( CI j = associations_.begin();
	     j != associations_.end(); j++) {
	  QDPIO::cerr << j->first << endl << flush;
	}
	
	return this->OnUnknownType(id);
      }
      else {
	return (i->second)();
      }
    }
        
    //! Create the object
    /*! 
     * \param id       object id
     * \param p1       first parameter arg to callback
     * \return returns pointer to the object
     */
    AbstractProduct* createObject(const IdentifierType& id, Parm1 p1)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) {
	
	typedef typename IdToProductMap::const_iterator CI;
	QDPIO::cerr << "Couldnt find key " << id << " in the map: " << endl;
	QDPIO::cerr << "Available Keys are : " << endl;
	for( CI j = associations_.begin();
	     j != associations_.end(); j++) {
	  QDPIO::cerr << j->first << endl << flush;
	}

	return this->OnUnknownType(id);
      }
      else {
	return (i->second)(p1);
      }
    }
        
    AbstractProduct* createObject(const IdentifierType& id, Parm1 p1, Parm2 p2)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) {

	typedef typename IdToProductMap::const_iterator CI;
	QDPIO::cerr << "Couldnt find key " << id << " in the map: " << endl;
	QDPIO::cerr << "Available Keys are : " << endl;
	for( CI j = associations_.begin();
	     j != associations_.end(); j++) {
	  QDPIO::cerr << j->first << endl << flush;
	}
	
	return this->OnUnknownType(id);
      }
      else {
	  return (i->second)(p1, p2);
      }
    }
    
    AbstractProduct* createObject(const IdentifierType& id, Parm1 p1, Parm2 p2, Parm3 p3)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) {
	
	typedef typename IdToProductMap::const_iterator CI;
	QDPIO::cerr << "Couldnt find key " << id << " in the map: " << endl;
	QDPIO::cerr << "Available Keys are : " << endl;
	for( CI j = associations_.begin();
	     j != associations_.end(); j++) {
	  QDPIO::cerr << j->first << endl << flush;
	}
	
	return this->OnUnknownType(id);
      }
      else {
	  return (i->second)(p1, p2, p3);
      }
    }
        
    AbstractProduct* createObject(const IdentifierType& id, Parm1 p1, Parm2 p2, Parm3 p3,
				  Parm4 p4)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) {
	
	typedef typename IdToProductMap::const_iterator CI;
	QDPIO::cerr << "Couldnt find key " << id << " in the map: " << endl;
	QDPIO::cerr << "Available Keys are : " << endl;
	for( CI j = associations_.begin();
	     j != associations_.end(); j++) {
	  QDPIO::cerr << j->first << endl << flush;
	}
	
	return this->OnUnknownType(id);
      }
      else {
	return (i->second)(p1, p2, p3, p4);
      }
    }
        
    AbstractProduct* createObject(const IdentifierType& id, Parm1 p1, Parm2 p2, Parm3 p3,
				  Parm4 p4, Parm5 p5)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) {
	
	typedef typename IdToProductMap::const_iterator CI;
	QDPIO::cerr << "Couldnt find key " << id << " in the map: " << endl;
	QDPIO::cerr << "Available Keys are : " << endl;
	for( CI j = associations_.begin();
	     j != associations_.end(); j++) {
	  QDPIO::cerr << j->first << endl << flush;
	}
	
	return this->OnUnknownType(id);
      }
      else {
	return (i->second)(p1, p2, p3, p4, p5);
      }
    }

    AbstractProduct* createObject(const IdentifierType& id, Parm1 p1, Parm2 p2, Parm3 p3,
				  Parm4 p4, Parm5 p5, Parm6 p6)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) {
	
	typedef typename IdToProductMap::const_iterator CI;
	QDPIO::cerr << "Couldnt find key " << id << " in the map: " << endl;
	QDPIO::cerr << "Available Keys are : " << endl;
	for( CI j = associations_.begin();
	     j != associations_.end(); j++) {
	  QDPIO::cerr << j->first << endl << flush;
	}
	
	return this->OnUnknownType(id);
      }
      else {
	return (i->second)(p1, p2, p3, p4, p5, p6);
      }
    }

    AbstractProduct* createObject(const IdentifierType& id, Parm1 p1, Parm2 p2, Parm3 p3,
				  Parm4 p4, Parm5 p5, Parm6 p6, Parm7 p7)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) {
	
	typedef typename IdToProductMap::const_iterator CI;
	QDPIO::cerr << "Couldnt find key " << id << " in the map: " << endl;
	QDPIO::cerr << "Available Keys are : " << endl;
	for( CI j = associations_.begin();
	     j != associations_.end(); j++) {
	  QDPIO::cerr << j->first << endl << flush;
	}
	
	return this->OnUnknownType(id);
      }
      else {
	return (i->second)(p1, p2, p3, p4, p5, p6, p7);
      }
    }

    AbstractProduct* createObject(const IdentifierType& id, Parm1 p1, Parm2 p2, Parm3 p3,
				  Parm4 p4, Parm5 p5, Parm6 p6, Parm7 p7, Parm8 p8)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) {
	
	typedef typename IdToProductMap::const_iterator CI;
	QDPIO::cerr << "Couldnt find key " << id << " in the map: " << endl;
	QDPIO::cerr << "Available Keys are : " << endl;
	for( CI j = associations_.begin();
	     j != associations_.end(); j++) {
	  QDPIO::cerr << j->first << endl << flush;
	}
	
	return this->OnUnknownType(id);
      }
      else {
	return (i->second)(p1, p2, p3, p4, p5, p6, p7, p8);
      }
    }

  
  private:
    typedef std::map<IdentifierType, ProductCreator> IdToProductMap;
    IdToProductMap associations_;
  };

}   // namespace Chroma

#endif
