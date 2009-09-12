#include "xml_help.h"
#include <DBFunc.h>
#include <string>
using namespace std;

namespace Chroma {
  namespace LaphEnv {


// *********************************************************

    // This returns the number of times that the tag "tagname"
    // is found in the descendents of the current context or the current context.
    // A **tag name** should be input, not an Xpath, since
    // a "./descendant-or-self::" is prepended to form an Xpath.

int xml_tag_count(XMLReader& xmlr, const string& tagname)
{
 return xmlr.count("./descendant-or-self::"+tagname);
}

 // ************************************************************

    //  This routine searches in the XML string "in_str" to find
    //  the tag with name "tag_name" and returns that complete
    //  XML element in "out_str".  If the tag is not found or
    //  multiple occurrences are found, this indicates a serious
    //  error and an abort is issued.  "callingClass" is the name
    //  of the class that called this function and is output 
    //  before an abort to help correct the situation.

void extract_xml_element(const string& in_str, const string& tag_name,
                         string& out_str, const string& callingClass)
{
 int start=in_str.find("<"+tag_name+">");
 int stop=in_str.find("</"+tag_name+">");
 if ((start==string::npos)||(stop==string::npos)||(stop<start)){
    QDPIO::cerr << "could not find tag "<<tag_name<<" in header string"
                << " in class "<<callingClass<<endl;
    QDP_abort(1);}
 out_str=in_str.substr(start,stop-start+tag_name.length()+3);
 if (in_str.find("<"+tag_name+">",start+1)!=string::npos){
    QDPIO::cerr << "multiple occurrences of tag "<<tag_name
                <<" in header string in class "<<callingClass<<endl;
    QDP_abort(1);}
}

// *********************************************************

   //  This routine checks that the current "top" node in "xmlr"
   //  has name "tagname".

bool checkXMLCurrentTagName(XMLReader& xmlr, const std::string& tagname)
{
 ostringstream path;
 path << "self::"<<tagname;
 return (xmlr.count(path.str())==1);
// stringstream oss;
// xmlr.print(oss);
// string xmlstr=oss.str();
// int start=xmlstr.find_first_of("<");
// if (start==string::npos) return false;
// int stop=xmlstr.find_first_of(">");
// if (stop==string::npos) return false;
// xmlstr=tidyString(xmlstr.substr(start+1,stop-start-1));
// return (xmlstr==tagname);
}

 // ************************************************************

  // removes tabs, newline, linefeed characters, then trims
  // leading and trailing blanks.

string tidyString(const string& str)   
{
 string tmp;
 for (int i=0;i<str.length();i++)
    if ((str[i]!='\n')&&(str[i]!='\t')&&(str[i]!='\r'))
       tmp.push_back(str[i]);
 int start=tmp.find_first_not_of(" ");
 if (start==string::npos) return "";
 int len=tmp.find_last_not_of(" ")-start+1;
 return tmp.substr(start,len);
}

// *************************************************************


bool xmlContentIsEqual(const string& doc1, const string& doc2, 
                       float float_rel_tol)
{
 bool result;
 if (Layout::primaryNode()){
    XMLContentComparer temp;
    result=temp.IsEqual(doc1,doc2,float_rel_tol);}
 Internal::broadcast(result);
 return result;  
}

bool xmlContentIsEqual(XMLReader& xmlr1, XMLReader& xmlr2, 
                       float float_rel_tol)
{
 stringstream oss1,oss2;
 xmlr1.print(oss1);
 xmlr2.print(oss2);
 xmlContentIsEqual(oss1.str(),oss2.str(),float_rel_tol);
}

bool headerMatch(const string& doc1, const string& doc2, 
                 float float_rel_tol)
{
 if (doc1==doc2) return true;
 return xmlContentIsEqual(doc1,doc2,float_rel_tol);
}

// ************************************************************

    //  tests if file having name "file_name" exists
    //  on the primary node or not; broadcasts results
    //  to all nodes

bool fileExists(const std::string& file_name)
{
 bool result;
 if (Layout::primaryNode()){
   result = FILEDB::fileExists(file_name);
 }
 Internal::broadcast(result);
 return result;
}


// ************************************************************



bool XML_NodePtr::operator<(const XML_NodePtr& rhs)
{
 return (m_ctl->node_cmp(m_ptr,rhs.m_ptr)<0);
}

bool XML_AttrPtr::operator<(const XML_AttrPtr& rhs)
{
 return (m_ctl->attr_cmp(m_ptr,rhs.m_ptr)<0);
}


bool XMLContentComparer::IsEqual(const string& doc1, const string& doc2, 
                                 float float_rel_tol)
{
 m_float_rel_tol=abs(float_rel_tol);
 if (m_float_rel_tol<1e-12) m_float_rel_tol=1e-12;

    // Parse XML document from a character string

 xmlDocPtr xml1 = xmlParseMemory(doc1.c_str(),doc1.size());
 if (xml1 == 0){
    QDPIO::cerr << "Error in XMLContentComparer::IsEqual"<<endl;
    QDPIO::cerr << "   ...could not parse first XML document"<<endl;
    xmlFreeDoc(xml1);
    return false;}

 xmlDocPtr xml2 = xmlParseMemory(doc2.c_str(),doc2.size());
 if (xml2 == 0){
    QDPIO::cerr << "Error in XMLContentComparer::IsEqual"<<endl;
    QDPIO::cerr << "   ...could not parse second XML document"<<endl;
    xmlFreeDoc(xml1);
    xmlFreeDoc(xml2);
    return false;}

    //Get the root element nodes
 xmlNodePtr root1 = xmlDocGetRootElement(xml1);
 xmlNodePtr root2 = xmlDocGetRootElement(xml2);

 int result=sibling_node_cmp(root1,root2);

 xmlFreeDoc(xml1);
 xmlFreeDoc(xml2);

 return (result==0);
}


  // Beginning with character at "start", moves the "character"
  // pointer to the first non-ignorable character (or NULL if
  // no such characters found).  "nchar" will contain the 
  // number of non-ignorable characters in the next token. 

void XMLContentComparer::get_next_token(const char*& start, int& nchar)
{
 char delimit[] = " \t\n\r";
 start += strspn(start, delimit); 
 if (*start=='\0'){ start=0; nchar=0; return;}
 const char *stop = strpbrk(start, delimit);
 if (stop==0) nchar=strlen(start);
 else{
    nchar=0; const char *tmp=start;
    while (tmp!=stop){ tmp++; nchar++;}}
}

   //  Generalization of strcmp but not using null characters
   //  as string ending; uses numbers of characters in na,nb
   //  Return 0 if na-character string in "a" equals nb-character
   //  string in "b".  Negative return if a-string precedes b-string,
   //  and >0 return if a-string comes after b-string.

int XMLContentComparer::token_strcmp(const char *a, int na, 
                                     const char *b, int nb)
{
 if (na==nb) return strncmp(a,b,na);
 else if (na>nb){
    int k=strncmp(a,b,nb);
    if (k!=0) return k;
    else return 1;}
 else{
    int k=strncmp(a,b,na);
    if (k!=0) return k;
    else return -1;}
}

   // Check if the "nchar"-character token starting at "a"
   // is a valid integer; if so, return integer value in "ivalue".
   // We assume the token has no blank characters.

bool XMLContentComparer::is_integer_token(const char *const a, int nchar, 
                                          int& ivalue)
{
 if (nchar==0) return false;
 const char* ap=a;
 int k=0;
 if ((*ap=='+')||(*ap=='-')){ap++; k++;}
 if (k==nchar) return false;
 while (k<nchar){
    if (!isdigit(*ap)) return false;
    ap++; k++;}
 stringstream oss; oss << "%" << nchar << "d";
 return sscanf(a,oss.str().c_str(),&ivalue);
}

   // Check if the "nchar"-character token starting at "a"
   // is a valid floating-point number; if so, return value in "fvalue".
   // We assume the token has no blank characters.

bool XMLContentComparer::is_float_token(const char *const a, int nchar,
                                        float& fvalue)
{
 if (nchar==0) return false;
 const char* ap=a;
 int k=0;
 if ((*ap=='+')||(*ap=='-')){ap++; k++;}
 if (k==nchar) return false;
 bool dotflag=false, eflag=false;
 while (k<(nchar-1)){
    if (*ap=='.'){
       if (dotflag) return false;
       dotflag=true;}
    else if ((*ap=='e')||(*ap=='E')){
       if (eflag) return false;
       dotflag=true;
       eflag=true;  ap++; k++;
       if ((*ap!='+')&&(*ap!='-')){ ap--; k--;}
       else if (k==nchar) return false;}
    else if (!isdigit(*ap)) return false;
    ap++; k++;
    }
 if (!((isdigit(*ap))||((*ap=='.')&&(!dotflag)))) return false; 
     // last character must be digit or . (if not dots yet)
 stringstream oss; oss << "%" << nchar << "g";
 return sscanf(a,oss.str().c_str(),&fvalue);
}

    // This routine compares the content in the token strings "a" and "b"
    // having "na" and "nb" number of characters.  First the token are
    // compared as integers, then as floats (to within "float_rel_tol"),
    // then finally as strings.  A zero value is returned if they are
    // equal, a number > 0 if a>b and a number <0 is a<b.

int XMLContentComparer::token_content_cmp(const char *a, int na, 
                                          const char *b, int nb)
{
 int ia,ib;
 float fa,fb;
 if (is_integer_token(a,na,ia)&&is_integer_token(b,nb,ib)){
    if (ia>ib) return 1;
    else if (ia==ib) return 0;
    else return -1;}
 else if (is_float_token(a,na,fa) && is_float_token(b,nb,fb)){
    float favg=0.5*abs(fa+fb);
    if ((fa-fb)>(favg*m_float_rel_tol)) return 1;
    else if ((fa-fb)<-(favg*m_float_rel_tol)) return -1;
    else return 0;}
 return token_strcmp(a,na,b,nb);
}


int XMLContentComparer::content_cmp(const char *a, const char *b)
{
 int result;
 const char *ap=a; int na;
 const char *bp=b; int nb;
 get_next_token(ap,na);
 get_next_token(bp,nb);
 while ((ap!=0)&&(bp!=0)){
    result=token_content_cmp(ap,na,bp,nb);
    if (result!=0) return result;
    ap+=na; get_next_token(ap,na);
    bp+=nb; get_next_token(bp,nb);
    }
 if ((ap==0)&&(bp==0)) return 0;
 if (ap!=0) return 1;
 else return -1;
}


int XMLContentComparer::attr_cmp(xmlAttrPtr a_attr, 
                                 xmlAttrPtr b_attr)
{
    // first compare by attribute name
 int result = strcmp((char*) a_attr->name, (char*) b_attr->name);
 if (result!=0) return result;

    // compare by value
 if  ((!xmlNodeIsText(a_attr->children))
    ||(!xmlNodeIsText(b_attr->children))) QDP_abort(1);  // quick check
 return strcmp((char*) a_attr->children->content,
               (char*) b_attr->children->content);
}


int XMLContentComparer::attrlist_cmp(xmlAttrPtr a_attr, 
                                     xmlAttrPtr b_attr)
{
 xmlAttrPtr cur_attr;

 list<XML_AttrPtr> alist;
 for (cur_attr=a_attr; cur_attr; cur_attr=cur_attr->next)
    alist.push_back(XML_AttrPtr(cur_attr,this));
 alist.sort();

 list<XML_AttrPtr> blist;
 for (cur_attr=b_attr; cur_attr; cur_attr=cur_attr->next)
    blist.push_back(XML_AttrPtr(cur_attr,this));
 blist.sort();

 int na=alist.size();
 int nb=blist.size();
 int n = (na>=nb)?nb:na;

 list<XML_AttrPtr>::const_iterator at=alist.begin();
 list<XML_AttrPtr>::const_iterator bt=blist.begin();

 for (int k=0;k<n;k++){
    int result=attr_cmp(at->m_ptr,bt->m_ptr);
    if (result!=0) return result;
    at++; bt++;}
 if ((at==alist.end())&&(bt==blist.end())) return 0;
 if (at==alist.end()) return -1;
 else return 1;
}

int XMLContentComparer::node_cmp(xmlNodePtr a_node,
                                 xmlNodePtr b_node)
{
    // compare node types first
 if (a_node->type < b_node->type) return -1;
 else if (a_node->type > b_node->type) return 1;

    // compare node names next
 int result=strcmp((char*) a_node->name, (char*) b_node->name);
 if (result!=0) return result;

    // compare by attributes thirdly
 result=attrlist_cmp(a_node->properties,b_node->properties);
 if (result!=0) return result;

    // if text node, compare by textual content
 if (xmlNodeIsText(a_node)){
    if (!xmlNodeIsText(b_node)) QDP_abort(1); // something went wrong
    return content_cmp((char*) a_node->content, (char*) b_node->content);}

    // now compare by children
 return sibling_node_cmp(a_node->children,b_node->children);
}


int XMLContentComparer::sibling_node_cmp(xmlNodePtr a_node,
                                         xmlNodePtr b_node)
{
 xmlNodePtr cur_node;

 list<XML_NodePtr> alist;
 for (cur_node=a_node; cur_node; cur_node=cur_node->next)
    alist.push_back(XML_NodePtr(cur_node,this));
 alist.sort();

 list<XML_NodePtr> blist;
 for (cur_node=b_node; cur_node; cur_node=cur_node->next)
    blist.push_back(XML_NodePtr(cur_node,this));
 blist.sort();

 int na=alist.size();
 int nb=blist.size();
 int n = (na>=nb)?nb:na;

 list<XML_NodePtr>::const_iterator at=alist.begin();
 list<XML_NodePtr>::const_iterator bt=blist.begin();

 for (int k=0;k<n;k++){
    int result=node_cmp(at->m_ptr,bt->m_ptr);
    if (result!=0) return result;
    at++; bt++;}
 if ((at==alist.end())&&(bt==blist.end())) return 0;
 if (at==alist.end()) return -1;
 else return 1;
}

// **********************************************************
  }
}
