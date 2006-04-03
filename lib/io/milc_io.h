// -*- C++ -*-
// $Id: milc_io.h,v 3.0 2006-04-03 04:58:55 edwards Exp $

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
  MILCGauge_t();          /*!< Default constructor */
  multi1d<int> nrow;      /*!< Lattice size */
  std::string  date;      /*!< ASCII date and time of file creation (64 bytes long) */
};


//! Source header read
void read(XMLReader& xml, const std::string& path, MILCGauge_t& header);

//! Source header writer
void write(XMLWriter& xml, const std::string& path, const MILCGauge_t& header);

}  // end namespace Chroma

#endif
