// $Id: cppacs_io.h,v 3.0 2006-04-03 04:58:55 edwards Exp $

/*! \file
 *  \brief CPPACS gauge format routines
 */

#ifndef __cppacs_io_h__
#define __cppacs_io_h__

#include "chromabase.h"

namespace Chroma {


//! CPPACS gauge field header
struct CPPACSGauge_t
{
  multi1d<int> nrow;      // Lattice size
  std::string  date;      // ASCII date and time of file creation 
};


//! Initialize header with default values
void CPPACSGaugeInit(CPPACSGauge_t& header);

//! Source header read
void read(XMLReader& xml, const std::string& path, CPPACSGauge_t& header);

//! Source header writer
void write(XMLWriter& xml, const std::string& path, const CPPACSGauge_t& header);

}  // end namespace Chroma

#endif
