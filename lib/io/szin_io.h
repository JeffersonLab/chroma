// $Id: szin_io.h,v 1.1 2003-10-08 02:30:26 edwards Exp $

/*! \file
 *  \brief  Routines associated with SZIN gauge field IO
 */

#ifndef __szin_io_h__
#define __szin_io_h__

//! Source header read
void read(XMLReader& xml, const string& path, SzinGauge_t& header);

//! Source header writer
void write(XMLWriter& xml, const string& path, const SzinGauge_t& header);

#endif
