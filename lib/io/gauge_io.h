// $Id: gauge_io.h,v 3.0 2006-04-03 04:58:55 edwards Exp $

/*! \file
 *  \brief Gauge reader/writers in QIO format
 */

#ifndef __gauge_io_h__
#define __gauge_io_h__

#include "chromabase.h"

namespace Chroma {

//! Read a gauge config in QIO format
/*!
 * \ingroup io
 *
 * \param file_xml     xml reader holding config info ( Modify )
 * \param record_xml   xml reader holding config info ( Modify )
 * \param u            gauge configuration ( Modify )
 * \param file         path ( Read )
 * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */    

void readGauge(XMLReader& file_xml, 
	       XMLReader& record_xml,
	       multi1d<LatticeColorMatrix>& u, 
	       const std::string& file, 
	       QDP_serialparallel_t serpar);


//! Write a gauge config in QIO format
/*!
 * \ingroup io
 *
 * \param file_xml    xml reader holding config info ( Modify )
 * \param record_xml  xml reader holding config info ( Modify )
 * \param u           gauge configuration ( Read )
 * \param cfg_file    path ( Read )
 * \param volfmt      either QDPIO_SINGLEFILE, QDPIO_MULTIFILE ( Read )
 * \param serpar      either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */    

void writeGauge(XMLBufferWriter& file_xml,
		XMLBufferWriter& record_xml, 
		const multi1d<LatticeColorMatrix>& u, 
		const std::string& file, 
		QDP_volfmt_t volfmt, 
		QDP_serialparallel_t serpar);

}  // end namespace Chroma

#endif
