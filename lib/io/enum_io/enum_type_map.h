#ifndef enum_type_map_h
#define enum_type_map_h

#include "chromabase.h"
#include <iostream>
#include <string>
#include <map>

using namespace std;
namespace Chroma { 

  template<typename EnumType>
    class EnumTypeMap {
    public:
    typedef typename std::map<std::string, EnumType> Str2Enum;
    typedef typename std::map<EnumType, std::string> Enum2Str;
    typedef typename std::map<string, EnumType>::iterator Str2EnumIterator;
    typedef typename std::map<EnumType, string>::iterator Enum2StrIterator;
    private:
    Enum2Str enum_map;  // Map one way (Enum->String)
    Str2Enum str_map;   // Map other way (String->Enum)
    EnumTypeMap(const EnumTypeMap<EnumType>& t) {}
    public:

    EnumTypeMap(){};

    // Function registers a string/enum pair
    bool registerPair(const string& s, const EnumType t) { 
      bool success = enum_map.insert(make_pair(t,s)).second;
      success &= str_map.insert(make_pair(s,t)).second;
      return success;
    }
  
    // Look up an enum based on a string 
    EnumType lookUpEnum(const string& s) {
      // Do the lookup
      Str2EnumIterator it = str_map.find(s);
      EnumType t;
      if( it != str_map.end() ) { 
	t= it->second;
      }
      else {      
	std::string s = "lookUpEnum failed on key: " + s;
	throw s;
      }   
      return t;
    }

    // Look up a string from an enum 
    string lookUpString(const EnumType& t) {

      // Do the lookup
      Enum2StrIterator it = enum_map.find(t);
      string s;
      if( it != enum_map.end() ) { 
	s= it->second;
      }
      else {      
	throw std::string("err");
      }   
      return s;
    }

    // "Reader"
    void read(XMLReader& xml_in, const string& path, EnumType& t)
    {
      string string_in;
      QDP::read(xml_in, path, string_in);
      t = lookUpEnum(string_in);
    }

    // Writer
    void write(XMLWriter& xml_out, const string& path, const EnumType& t) {
      string string_out = lookUpString(t);
      QDP::write(xml_out, path, string_out);
    }
  };	      
  

};


#endif 
