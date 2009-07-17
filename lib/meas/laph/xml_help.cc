#include "xml_help.h"
#include <string>
using namespace std;

namespace Chroma {
  namespace LaphEnv {

 // ************************************************************

   //  This routine checks that the current "top" node in "xmlr"
   //  has name "tagname".

bool checkXMLCurrentTagName(XMLReader& xmlr, const std::string& tagname)
{
 stringstream oss;
 xmlr.print(oss);
 string xmlstr=oss.str();
 int start=xmlstr.find_first_of("<");
 if (start==string::npos) return false;
 int stop=xmlstr.find_first_of(">");
 if (stop==string::npos) return false;
 xmlstr=tidyString(xmlstr.substr(start+1,stop-start-1));
 return (xmlstr==tagname);
}

 // ************************************************************

  // removes tabs and newline characters, then trims
  // leading and trailing blanks.

string tidyString(const string& str)   
{
 string tmp;
 for (int i=0;i<str.length();i++)
    if ((str[i]!='\n')&&(str[i]!='\t'))
       tmp.push_back(str[i]);
 int start=tmp.find_first_not_of(" ");
 if (start==string::npos) return "";
 int len=tmp.find_last_not_of(" ")-start+1;
 return tmp.substr(start,len);
}

// **********************************************************
  }
}
