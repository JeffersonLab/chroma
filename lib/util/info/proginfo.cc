//  $Id: proginfo.cc,v 1.1 2004-01-29 20:38:09 edwards Exp $
/*! \file
 *  \brief Print out basic info about this program
 */

#include "chromabase.h"
#include "util/geom/printgeom.h"
#include "util/geom/proginfo.h"

using namespace QDP;

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
  START_CODE("proginfo");

  push(xml,"ProgramInfo");

  write(xml,"run_date",date);
  printgeom(xml);

  pop(xml);

  END_CODE("proginfo");
}
