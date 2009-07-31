#ifndef XML_HELP_H
#define XML_HELP_H

#include "qdp.h"
#include "chromabase.h"
#include <libxml/tree.h>

namespace Chroma {
  namespace LaphEnv {

// *********************************************************

   //  This routine checks that the current "top" node in "xmlr"
   //  has name "tagname".

bool checkXMLCurrentTagName(XMLReader& xmlr, const std::string& tagname);


// *********************************************************


  // Removes tabs and newline characters, then trims
  // leading and trailing blanks.

std::string tidyString(const std::string& str);    


// *********************************************************


  // Compares the XML content in "doc1" and "doc2" and returns
  // "true" if they are the same.  In doing the comparison, the
  // order of sibling nodes is irrelevant, and textual content
  // is compared token by token.  Integer tokens are compared as
  // integers, and floating-point tokens are compared as 
  // floats.  If the difference between floats is less than
  // "float_rel_tol", they are considered the same.


bool xmlContentIsEqual(const string& doc1, const string& doc2, 
                       float float_rel_tol = 1e-6);


  // Same as above, but applied to XMLReaders "xmlr1" and "xmlr2".

bool xmlContentIsEqual(XMLReader& xmlr1, XMLReader& xmlr2, 
                       float float_rel_tol = 1e-6);


// *********************************************************

    //  tests if file having name "file_name" exists
    //  on the primary node or not; broadcasts results
    //  to all nodes

bool fileExists(const std::string& file_name);



// *********************************************************


  // declarations below are NOT available to the end user;
  // these classes are used by xmlContentIsEqual

class XMLContentComparer;  // forward declaration


class XML_NodePtr {

    xmlNodePtr m_ptr;
    XMLContentComparer *m_ctl;

    XML_NodePtr() {}
    XML_NodePtr(xmlNodePtr ptr, XMLContentComparer *ctl) 
            : m_ptr(ptr), m_ctl(ctl) {}

 public:

    XML_NodePtr(const XML_NodePtr& in) 
            : m_ptr(in.m_ptr), m_ctl(in.m_ctl)  {}

    ~XML_NodePtr() {}

    friend class XMLContentComparer;

    bool operator<(const XML_NodePtr& rhs);

};

class XML_AttrPtr {

    xmlAttrPtr m_ptr;
    XMLContentComparer *m_ctl;

    XML_AttrPtr() {}
    XML_AttrPtr(xmlAttrPtr ptr, XMLContentComparer *ctl) 
            : m_ptr(ptr), m_ctl(ctl) {}

 public:

    XML_AttrPtr(const XML_AttrPtr& in) 
            : m_ptr(in.m_ptr), m_ctl(in.m_ctl)  {}

    ~XML_AttrPtr() {}

    friend class XMLContentComparer;

    bool operator<(const XML_AttrPtr& rhs);

};

class XMLContentComparer {

    float m_float_rel_tol;

       // for use only by xmlContentIsEqual
    XMLContentComparer() {}
    XMLContentComparer(const XMLContentComparer& in);
    XMLContentComparer& operator=(const XMLContentComparer& in);


 public:

    ~XMLContentComparer() {}

 private:

    friend bool xmlContentIsEqual(const string& doc1, const string& doc2, 
                                  float float_rel_tol);
    friend class XML_NodePtr;
    friend class XML_AttrPtr;

    bool IsEqual(const string& doc1, const string& doc2, 
                 float float_rel_tol);

    void get_next_token(const char*& start, int& nchar);

    int token_strcmp(const char *a, int na, const char *b, int nb);

    bool is_integer_token(const char *const a, int nchar, int& ivalue);

    bool is_float_token(const char *const a, int nchar, float& fvalue);

    int token_content_cmp(const char *a, int na, const char *b, int nb);

    int content_cmp(const char *a, const char *b);

    int attr_cmp(xmlAttrPtr a_attr, xmlAttrPtr b_attr);

    int attrlist_cmp(xmlAttrPtr a_attr, xmlAttrPtr b_attr);

    int sibling_node_cmp(xmlNodePtr a_node, xmlNodePtr b_node);

    int node_cmp(xmlNodePtr a_node, xmlNodePtr b_node);


};



// *********************************************************
  }
}
#endif
