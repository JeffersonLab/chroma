//  $Id: proginfo.cc,v 1.2 2004-01-29 20:46:29 edwards Exp $
/*! \file
 *  \brief Print out basic info about this program
 */

#include "chromabase.h"
#include "util/info/printgeom.h"
#include "util/info/proginfo.h"

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

//  write(xml,"run_date",date);
  printgeom(xml);

  pop(xml);

  END_CODE("proginfo");
}
