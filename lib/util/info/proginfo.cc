//  $Id: proginfo.cc,v 1.6 2005-01-31 17:56:15 edwards Exp $
/*! \file
 *  \brief Print out basic info about this program
 */

#include "chromabase.h"
#include "util/info/printgeom.h"
#include "util/info/proginfo.h"
#include <ctime>


namespace Chroma {

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
    write(xml, CHROMA_PACKAGE_NAME, CHROMA_PACKAGE_VERSION);
    write(xml, "qdp", QDP_PACKAGE_VERSION);
    pop(xml);

    std::time_t now;
    if(std::time(&now)==-1)
	QDPIO::cerr<<"proginfo: Cannot get the time.\n";
    std::tm *tp = std::localtime(&now);

    char date[64];
    std::strftime(date, 63, "%d %b %y %X %Z", tp);
    write(xml,"run_date", date);

    printgeom(xml);

    pop(xml);

    END_CODE();
}

}  // end namespace Chroma
