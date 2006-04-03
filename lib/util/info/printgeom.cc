//  $Id: printgeom.cc,v 3.0 2006-04-03 04:59:13 edwards Exp $
/*! \file
 *  \brief Print out machine geometry and problem size info
 */

#include "chromabase.h"
#include "util/info/printgeom.h"

namespace Chroma {

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
  START_CODE();

  push(xml,"Setgeom");
  write(xml,"latt_size",Layout::lattSize());
  write(xml,"logical_size",Layout::logicalSize());
  write(xml,"subgrid_size",Layout::subgridLattSize());
  write(xml,"total_volume",Layout::vol());
  write(xml,"subgrid_volume",Layout::sitesOnNode());
  pop(xml);

  END_CODE();
}

}  // end namespace Chroma
