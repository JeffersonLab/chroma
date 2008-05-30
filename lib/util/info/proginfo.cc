//  $Id: proginfo.cc,v 3.2 2008-05-30 18:53:11 edwards Exp $
/*! \file
 *  \brief Print out basic info about this program
 */

#include "chroma_config.h"
#include "chromabase.h"
#include "util/info/printgeom.h"
#include "util/info/proginfo.h"
#include <ctime>


namespace Chroma  
{

  //! Print out basic information about this program
  /*!
   * \ingroup info
   *
   *
   *  \param xml          The XML stream to which to write the information.
   */

  void proginfo(XMLWriter& xml)
  {
    START_CODE();

    push(xml,"ProgramInfo");

    push(xml,"code_version");
    write(xml, "chroma", CHROMA_PACKAGE_VERSION);
    write(xml, "qdp", QDP_PACKAGE_VERSION);
    write(xml, "basePrecision", BASE_PRECISION);
    pop(xml);

    std::time_t now;

    if(std::time(&now)==-1)
    {
      QDPIO::cerr<<"proginfo: Cannot get the time.\n";
      return;
    }
    std::tm *tp = std::localtime(&now);

    char date[64];
    std::strftime(date, 63, "%d %b %y %X %Z", tp);
    write(xml,"run_date", date);

    printgeom(xml);
    pop(xml);

    END_CODE();
  }

}  // end namespace Chroma
