//  $Id: proginfo.cc,v 1.4 2005-01-14 18:42:38 edwards Exp $
/*! \file
 *  \brief Print out basic info about this program
 */

#include "chromabase.h"
#include "util/info/printgeom.h"
#include "util/info/proginfo.h"

namespace Chroma {

//! Print out basic info about this program
/*!
 * \ingroup info
 *
 * Arguments:
 *
 *  \param xml          The xml stream to write the info
 */

void proginfo(XMLWriter& xml)
{
  START_CODE();

  push(xml,"ProgramInfo");

//  write(xml,"run_date",date);
  printgeom(xml);

  pop(xml);

  END_CODE();
}

}  // end namespace Chroma
