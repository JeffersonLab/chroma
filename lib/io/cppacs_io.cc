// $Id: cppacs_io.cc,v 1.1 2005-01-13 16:27:58 mcneile Exp $

/*! \file
 *  \brief CPPACS gauge format routines
 */

#include "chromabase.h"
#include "io/cppacs_io.h"

using namespace QDP;

//! Initialize header with default values
void CPPACSGaugeInit(CPPACSGauge_t& header)
{
  header.nrow = Layout::lattSize();
  header.date = "put in some date here";
}



//! Source header read
void read(XMLReader& xml, const string& path, CPPACSGauge_t& header)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "date", header.date);
  read(paramtop, "nrow", header.nrow);
}


//! Source header writer
void write(XMLWriter& xml, const string& path, const CPPACSGauge_t& header)
{
  push(xml, path);

  write(xml, "date", header.date);
  write(xml, "nrow", header.nrow);

  pop(xml);
}

