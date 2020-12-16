// -*- C++ -*-
/*! \file
 * \brief Enum std::map
 *
 * Enum std::map
 */

#ifndef enum_type_map_h
#define enum_type_map_h

#include "chromabase.h"
#include <iostream>
#include <map>

namespace Chroma 
{ 

  //! Main enum std::map holder
  /*! 
   * \ingroup io
   *
   * Enums used throughout the code and their initialization
   * via std::maps. This is used for IO.
   */
  template<typename EnumType>
  class EnumTypeMap 
  {
  public:
    typedef typename std::map<std::string, EnumType> Str2Enum;
    typedef typename std::map<EnumType, std::string> Enum2Str;
    typedef typename std::map<std::string, EnumType>::iterator Str2EnumIterator;
    typedef typename std::map<EnumType, std::string>::iterator Enum2StrIterator;
  private:
    Enum2Str enum_map;  // Map one way (Enum->String)
    Str2Enum str_map;   // Map other way (String->Enum)
    EnumTypeMap(const EnumTypeMap<EnumType>& t) {}

    struct EnumLookupException {
      EnumLookupException(const std::string& s) : val_str(s) {};
      const std::string val_str;
    };
    
    struct StringLookupException { 
      StringLookupException(const EnumType& t) : val_enum(t) {};
      const EnumType val_enum;
    };

    void dumpMapStrings(const std::string& typeIDString) {
      QDPIO::cout << "Allowed enum " << typeIDString << " std::string-value pairs are: " << std::endl;
      QDPIO::cout << "==================== " << std::endl;
      Enum2StrIterator it; 
      for( it=enum_map.begin(); it!=enum_map.end(); it++) {
	QDPIO::cout << "  String: " << it->second << ", Int value of enum: " << it->first << std::endl;
      }
    }
     
  public:

    EnumTypeMap(){};

    //! Function registers a std::string/enum pair
    /*!
     * Enum-String pairs are initialised when the relevant
     * env's "registered" boolean is set, and this may well be 
     * before QDP is initialised. So it is not safe to 
     * call QDPIO::cerr in this pair registration function....
     * Printing to std::cerr works tho. So I add it. Message
     * should only be generated if registration fails which
     * is most likely to be due to a duplicate key. This is 
     * a useful check that I haven't left duplicates in the 
     * registerAll() functions.
     */
    bool registerPair(const std::string& s, const EnumType t) { 
      bool success;
      success = enum_map.insert(make_pair(t,s)).second;
      if( ! success ) { 
	std::cerr << "Enum-String Insertion Failure. Enum="<<t<< " String=" << s << std::endl;
	std::cerr << "Most Likely pair already inserted. Check relevant files for duplicates" << std::endl << std::flush ;
	
	exit(1);
      }
      
      success &= str_map.insert(make_pair(s,t)).second;
      if( ! success ) { 
	std::cerr << "String-Enum Insertion Failure. Enum="<<t<< " String=" << s << std::endl;
	std::cerr << "Most Likely pair already inserted. Check relevant files for duplicates" << std::endl << std::flush ;
	
	exit(1);
      }
      
      return true;
    }
   

    //! Look up an enum based on a std::string 
    EnumType lookUpEnum(const std::string& s) {

      // Do the lookup
      Str2EnumIterator it = str_map.find(s);
      EnumType t;
      if( it != str_map.end() ) { 
	t= it->second;
      }
      else {     
	EnumLookupException e(s);
	throw e;
      }   
      return t;
    }
   

    //! Look up a std::string from an enum 
    std::string lookUpString(const EnumType& t) {
      // Do the lookup
      Enum2StrIterator it = enum_map.find(t);
      std::string s;
      if( it != enum_map.end() ) { 
	s= it->second;
      }
      else {      
	StringLookupException e(t);
	throw e;	
      }   
      return s;
    }
 

    //! "Reader"
    void read(const std::string& typeIDString, 
	      XMLReader& xml_in, 
	      const std::string& path, 
	      EnumType& t) {
      // Try to read a std::string for the enum
      std::string string_in;
      try { 
	QDP::read(xml_in, path, string_in);
      }
      catch( const std::string& xml_e ) { 
	QDPIO::cerr << "Caught Exception parsing XML: " << xml_e << std::endl;
	QDP_abort(1);
      }

      // Try to look up the enum for the std::string
      try { 
	t = lookUpEnum(string_in);
      }
      catch( EnumLookupException& e ) {
	QDPIO::cerr << "No enum " << typeIDString << " registered for std::string " << e.val_str << " while parsing XML Query " << path << std::endl;
	dumpMapStrings(typeIDString);
	QDP_abort(1);
      }
    }
   

    //! Writer
    void write(const std::string& typeIDString,
	       XMLWriter& xml_out, 
	       const std::string& path, 
	       const EnumType& t)  {
      std::string string_out;
      // Try to look up the std::string for the enum
      try { 
	string_out = lookUpString(t);
      }
      catch( StringLookupException& e) { 
	QDPIO::cerr << "No std::string found for enum " << typeIDString << " with value " << t << " while attempting to write XML " << std::endl;
	dumpMapStrings(typeIDString);
	QDP_abort(1);
      }

      // Try writing the std::string
      try { 
	QDP::write(xml_out, path, string_out);
      }
      catch( const std::string& xml_e ) { 
	QDPIO::cerr << "Caught Exception writing XML: " << xml_e << std::endl;
	QDP_abort(1);
      }
    }
    
    //! "Reader"
    void read(const std::string& typeIDString, 
	      BinaryReader& bin_in, 
	      EnumType& t) {
      // Try to read a std::string for the enum
      std::string string_in;
      try { 
	QDP::readDesc(bin_in, string_in);
      }
      catch( const std::string& bin_e ) { 
	QDPIO::cerr << "Caught Exception parsing BIN: " << bin_e << std::endl;
	QDP_abort(1);
      }

      // Try to look up the enum for the std::string
      try { 
	t = lookUpEnum(string_in);
      }
      catch( EnumLookupException& e ) {
	QDPIO::cerr << "No enum " << typeIDString << " registered for std::string " << e.val_str << " while parsing Binary lookup" << std::endl;
	dumpMapStrings(typeIDString);
	QDP_abort(1);
      }
    }
   

    //! Writer
    void write(const std::string& typeIDString,
	       BinaryWriter& bin_out, 
	       const EnumType& t)  {
      std::string string_out;
      // Try to look up the std::string for the enum
      try { 
	string_out = lookUpString(t);
      }
      catch( StringLookupException& e) { 
	QDPIO::cerr << "No std::string found for enum " << typeIDString << " with value " << t << " while attempting to write BIN " << std::endl;
	dumpMapStrings(typeIDString);
	QDP_abort(1);
      }

      // Try writing the std::string
      try { 
	QDP::writeDesc(bin_out, string_out);
      }
      catch( const std::string& bin_e ) { 
	QDPIO::cerr << "Caught Exception writing BIN: " << bin_e << std::endl;
	QDP_abort(1);
      }
    }
    
  };	      
  

}


#endif 
