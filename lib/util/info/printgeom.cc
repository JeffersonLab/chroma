//  $Id: printgeom.cc,v 1.2 2004-01-29 20:46:29 edwards Exp $
/*! \file
 *  \brief Print out machine geometry and problem size info
 */

#include "chromabase.h"
#include "util/info/printgeom.h"

using namespace QDP;

//! Print out machine geometry and problem size info
/*!
 * \ingroup info
 *
 * Arguments:
 *
 *  \param xml          The xml stream to write the info
 */

void printgeom(XMLWriter& xml)
{
  START_CODE("printgeom");

  push(xml,"Setgeom");
  write(xml,"latt_size",Layout::lattSize());
  write(xml,"logical_size",Layout::logicalSize());
  write(xml,"subgrid_size",Layout::subgridLattSize());
  pop(xml);

  END_CODE("printgeom");
}
