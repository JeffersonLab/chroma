// $Id: milc_io.h,v 1.2 2005-01-14 20:13:06 edwards Exp $

/*! \file
 *  \brief MILC gauge format routines
 */

#ifndef __milc_io_h__
#define __milc_io_h__

#include "chromabase.h"

namespace Chroma {


//! MILC gauge field header
struct MILCGauge_t
{
  multi1d<int> nrow;      // Lattice size
  std::string  date;      // ASCII date and time of file creation (64 bytes long)
};


//! Initialize header with default values
void MILCGaugeInit(MILCGauge_t& header);

//! Source header read
void read(XMLReader& xml, const std::string& path, MILCGauge_t& header);

//! Source header writer
void write(XMLWriter& xml, const std::string& path, const MILCGauge_t& header);

}  // end namespace Chroma

#endif
