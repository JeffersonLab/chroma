
/*! \file
 *  \brief CPPACS gauge format routines
 */

#include "chromabase.h"
#include "io/cppacs_io.h"

namespace Chroma {


//! Initialize header with default values
void CPPACSGaugeInit(CPPACSGauge_t& header)
{
  header.nrow = Layout::lattSize();
  header.date = "put in some date here";
}



//! Source header read
void read(XMLReader& xml, const std::string& path, CPPACSGauge_t& header)
{
  XMLReader paramtop(xml, path);

  read(paramtop, "date", header.date);
  read(paramtop, "nrow", header.nrow);
}


//! Source header writer
void write(XMLWriter& xml, const std::string& path, const CPPACSGauge_t& header)
{
  push(xml, path);

  write(xml, "date", header.date);
  write(xml, "nrow", header.nrow);

  pop(xml);
}

}  // end namespace Chroma
