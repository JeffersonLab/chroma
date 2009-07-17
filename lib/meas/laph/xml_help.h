#ifndef XML_HELP_H
#define XML_HELP_H

#include "qdp.h"
#include "chromabase.h"

namespace Chroma {
  namespace LaphEnv {

// *********************************************************

   //  This routine checks that the current "top" node in "xmlr"
   //  has name "tagname".

bool checkXMLCurrentTagName(XMLReader& xmlr, const std::string& tagname);


  // Removes tabs and newline characters, then trims
  // leading and trailing blanks.

std::string tidyString(const std::string& str);    


// *********************************************************
  }
}
#endif
