// $Id: gauge_io.cc,v 1.1 2004-04-05 19:48:19 edwards Exp $
/*! \file
 * \brief Routines associated with Chroma propagator gauge IO
 */

#include "chromabase.h"
#include "io/gauge_io.h"

using std::string;

// Read a Chroma propagator
/*
 * \param file_xml     xml reader holding config info ( Modify )
 * \param record_xml   xml reader holding config info ( Modify )
 * \param u            gauge configuration ( Modify )
 * \param file         path ( Read )
 * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */    
void readGauge(XMLReader& file_xml,
	       XMLReader& record_xml, 
	       multi1d<LatticeColorMatrix>& u, 
	       const string& file, 
	       QDP_serialparallel_t serpar)
{
  QDPFileReader to(file_xml,file,serpar);
  read(to,record_xml,u);
  close(to);
}


// Write a Gauge field in QIO format
/*
 * \param file_xml    xml reader holding config info ( Modify )
 * \param record_xml  xml reader holding config info ( Modify )
 * \param u           gauge configuration ( Modify )
 * \param file        path ( Read )
 * \param volfmt      either QDPIO_SINGLEFILE, QDPIO_MULTIFILE ( Read )
 * \param serpar      either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */    
void writeGauge(XMLBufferWriter& file_xml,
		XMLBufferWriter& record_xml, 
		const multi1d<LatticeColorMatrix>& u, 
		const string& file, 
		QDP_volfmt_t volfmt, 
		QDP_serialparallel_t serpar)
{
  QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
  write(to,record_xml,u);
  close(to);
}




