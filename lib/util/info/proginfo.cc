//  $Id: proginfo.cc,v 1.3 2004-07-28 02:38:06 edwards Exp $
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
  START_CODE();

  push(xml,"ProgramInfo");

//  write(xml,"run_date",date);
  printgeom(xml);

  pop(xml);

  END_CODE();
}
