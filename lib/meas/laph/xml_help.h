#ifndef XML_HELP_H
#define XML_HELP_H

#include "qdp.h"
#include "chromabase.h"
#include <libxml/tree.h>

namespace Chroma {
  namespace LaphEnv {

// *********************************************************

    // This routine calls an xml read inside a try block
    // and outputs an informative message if the read fails.
    // The added parameter is a string which should be the
    // name of the class which called this function.
    // Note that a **tag name** should be used. A "./descendant-or-self::"
    // is prepended to form an appropriate Xpath.

template <typename T>
void xmlread(XMLReader& xmlr, const string& tagname, T& val,
             const string& callingClass)
{
 try{
    read(xmlr,"./descendant-or-self::"+tagname,val);}
 catch(const string& err_msg){
    QDPIO::cerr << "Invalid read of "<<tagname<<" in "<<callingClass<<endl;
    QDP_abort(1);}
}

// *********************************************************

    // This returns the number of times that the tag "tagname"
    // is found in the descendents of the current context.
    // A **tag name** should be input, not an Xpath, since
    // a "./descendant-or-self::" is prepended to form an Xpath.

int xml_tag_count(XMLReader& xmlr, const string& tagname);


// *********************************************************

template <typename T>
void assertEqual(const T& obj1, const T& obj2, const string& callingClass)
{
 try{
    obj1.checkEqual(obj2);}
 catch(const string& err_msg){
    QDPIO::cerr << err_msg <<" in "<<callingClass<<endl;
    QDP_abort(1);}
}
 
// *********************************************************

    //  This routine searches in the XML string "in_str" to find
    //  the tag with name "tag_name" and returns that complete
    //  XML element in "out_str".  If the tag is not found or
    //  multiple occurrences are found, this indicates a serious
    //  error and an abort is issued.  "callingClass" is the name
    //  of the class that called this function and is output 
    //  before an abort to help correct the situation.

void extract_xml_element(const string& in_str, const string& tag_name,
                         string& out_str, const string& callingClass);

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

  // Same as above, but first tries a straight string comparison.
  // If strings are the same, returns true; otherwise, an
  // XML content comparison is made.

bool headerMatch(const string& doc1, const string& doc2, 
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
