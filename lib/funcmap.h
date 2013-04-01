// -*- C++ -*-
// $Id: funcmap.h,v 3.0 2006-04-03 04:58:44 edwards Exp $
/*! @file
 * @brief Map of function calls
 */

#ifndef __funcmap_h__
#define __funcmap_h__

#include "chromabase.h"

#include "typeinfo.h"
#include "typelist.h"
#include <map>
#include <exception>


namespace Chroma
{

  ////////////////////////////////////////////////////////////////////////////////
  // class template DefaultFunctionMapError
  // Manages the "Unknown Type" error in an object factory
  ////////////////////////////////////////////////////////////////////////////////

  template <typename IdentifierType, class AbstractProduct>
  struct DefaultFunctionMapError
  {
    struct Exception : public std::exception
    {
      const char* what() const throw() { return "Unknown Type"; }
    };
        
    static AbstractProduct OnUnknownType(const IdentifierType& id)
    {
      throw Exception();
    }
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class template StringFunctionMapError
  // Manages the "Unknown Type" error in an object factory
  ////////////////////////////////////////////////////////////////////////////////

  template <typename IdentifierType, class AbstractProduct>
  struct StringFunctionMapError
  {
    static AbstractProduct OnUnknownType(const IdentifierType& id)
    {
      std::ostringstream os;
      os << "FunctionMap error: unknown identifier: id = " << id << endl;
      throw os.str();
    }
  };

  ////////////////////////////////////////////////////////////////////////////////
  // class template NullFunctionMapError
  // Manages the "Unknown Type" error in an object factory
  ////////////////////////////////////////////////////////////////////////////////

  template <typename IdentifierType, class AbstractProduct>
  struct NullFunctionMapError
  {
    static AbstractProduct OnUnknownType(const IdentifierType&)
    {
      return zero;
    }
  };


  ////////////////////////////////////////////////////////////////////////////////
  // class template DefaultDisambiguator
  ////////////////////////////////////////////////////////////////////////////////

  struct DefaultDisambiguator {};


  //! Object factory class
  /*! @ingroup actions
   *
   * Supports abstract creation of objects
   */
  template<typename Disambiguator,
           class AbstractProduct, 
	   typename IdentifierType,
           class TList = NullType,
           typename ProductCall = AbstractProduct (*)(),
           template<typename, class> class FunctionMapErrorPolicy = DefaultFunctionMapError>
  class FunctionMap : public FunctionMapErrorPolicy<IdentifierType, AbstractProduct>
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

    // Member functions
    //! Register the function
    /*! 
     * \param id       function id
     * \param call     the calling function
     * \return returns true if registration successful
     */
    bool registerFunction(const IdentifierType& id, ProductCall call)
    {
      return associations_.insert(
				  typename IdToProductMap::value_type(id, call)).second;
    }
  
    //! Unregister the function
    /*! 
     * \param id       function id
     * \return returns true if function was registered before
     */
    bool unregisterFunction(const IdentifierType& id)
    {
      return associations_.erase(id) == 1;
    }

    //! Call the function
    /*! 
     * \param id       function id
     * \return returns result of the function call
     */
    AbstractProduct callFunction(const IdentifierType& id)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) 
      {
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
        
    //! Call the function
    /*! 
     * \param id       function id
     * \param p1       first parameter arg to function
     * \return returns result of the function call
     */
    AbstractProduct callFunction(const IdentifierType& id, Parm1 p1)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) 
      {
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
        
    //! Call the function
    /*! 
     * \param id       function id
     * \param p1       1st parameter arg to function
     * \param p2       2nd parameter arg to function
     * \return returns result of the function call
     */
    AbstractProduct callFunction(const IdentifierType& id, Parm1 p1, Parm2 p2)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) 
      {
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
    
    //! Call the function
    /*! 
     * \param id       function id
     * \param p1       1st parameter arg to function
     * \param p2       2nd parameter arg to function
     * \param p3       3rd parameter arg to function
     * \return returns result of the function call
     */
    AbstractProduct callFunction(const IdentifierType& id, Parm1 p1, Parm2 p2, Parm3 p3)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) 
      {
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
        
    //! Call the function
    /*! 
     * \param id       function id
     * \param p1       1st parameter arg to function
     * \param p2       2nd parameter arg to function
     * \param p3       3rd parameter arg to function
     * \param p4       4th parameter arg to function
     * \return returns result of the function call
     */
    AbstractProduct callFunction(const IdentifierType& id, Parm1 p1, Parm2 p2, Parm3 p3,
				 Parm4 p4)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) 
      {
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
        
    //! Call the function
    /*! 
     * \param id       function id
     * \param p1       1st parameter arg to function
     * \param p2       2nd parameter arg to function
     * \param p3       3rd parameter arg to function
     * \param p4       4th parameter arg to function
     * \param p5       5th parameter arg to function
     * \return returns result of the function call
     */
    AbstractProduct callFunction(const IdentifierType& id, Parm1 p1, Parm2 p2, Parm3 p3,
				 Parm4 p4, Parm5 p5)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) 
      {
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

    //! Call the function
    /*! 
     * \param id       function id
     * \param p1       1st parameter arg to function
     * \param p2       2nd parameter arg to function
     * \param p3       3rd parameter arg to function
     * \param p4       4th parameter arg to function
     * \param p5       5th parameter arg to function
     * \param p6       6th parameter arg to function
     * \return returns result of the function call
     */
    AbstractProduct callFunction(const IdentifierType& id, Parm1 p1, Parm2 p2, Parm3 p3,
				 Parm4 p4, Parm5 p5, Parm6 p6)
    {
      typename IdToProductMap::const_iterator i = associations_.find(id);
      if (i == associations_.end()) 
      {
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

  
  private:
    typedef std::map<IdentifierType, ProductCall> IdToProductMap;
    IdToProductMap associations_;
  };

}   // namespace Chroma

#endif
