// $Id: readmilc.cc,v 1.2 2003-10-09 16:07:57 edwards Exp $

/*! \file
 *  \brief Read a MILC gauge configuration written in the 1997 format
 */

#include "chromabase.h"
#include "io/milc_io.h"
#include "io/readmilc.h"
#include "qdp_util.h"    // from QDP

#include <string>
using std::string;

using namespace QDP;

//! Read a MILC configuration file
/*!
 * \ingroup io
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readMILC(MILCGauge_t& header, multi1d<LatticeColorMatrix>& u, const string& cfg_file)
{
  START_CODE("readMILC");

  MILCGaugeInit(header);  // initialize the header with defaults

  BinaryReader cfg_in(cfg_file); // for now, cfg_io_location not used

  int magic_number;
  read(cfg_in, magic_number);

  if( magic_number != 20103)
    QDP_error_exit("readMILC: unexpected byte order of file");

  // Check lattice size
  header.nrow.resize(Nd);
  read(cfg_in, header.nrow);

  for(int j = 0; j < Nd; ++j)
    if ( header.nrow[j] != Layout::lattSize()[j] )
      QDP_error_exit("readMILC: unexpected lattice size: header.nrow[%d]=%d",
                     j,header.nrow[j]);

  // Time stamp
  char date_tmp[65];
  cfg_in.readArray(date_tmp, 64, 1);
  date_tmp[64] = '\0';
  header.date = date_tmp;

  // Site order - only support non-sitelist format
  int order;
  read(cfg_in, order);
  if( order != 0)
    QDP_error_exit("readMILC: only support non-sitelist format");

  /*
   * Read away...
   */
  
  // MILC format has the directions inside the sites
  multi1d<ColorMatrix> u_old(Nd);

  for(int site=0; site < Layout::vol(); ++site)
  {
    multi1d<int> coord = crtesn(site, Layout::lattSize()); // The coordinate
      
    // Read in Nd SU(3) matrices. 
    // NOTE: the su3_matrix layout should be the same as in QDP
    read(cfg_in, u_old); 

    for(int mu=0; mu < Nd; ++mu)
      pokeSite(u[mu], u_old[mu], coord); // Put it into the correct place
  }

  // Go ahead an read checksums, but will not use for now
  unsigned int sum29, sum31;
  read(cfg_in, sum29);
  read(cfg_in, sum31);

  cfg_in.close();

  END_CODE("readMILC");
}



//! Read a MILC configuration file
/*!
 * \ingroup io
 *
 * \param xml        xml reader holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readMILC(XMLReader& xml, multi1d<LatticeColorMatrix>& u, const string& cfg_file)
{
  START_CODE("readMILC");

  MILCGauge_t header;

  // Read the config and its binary header
  readMILC(header, u, cfg_file);

  // Now, set up the XML header. Do this by first making a buffer
  // writer that is then used to make the reader
  XMLBufferWriter  xml_buf;
  write(xml_buf, "MILC", header);

  try 
  {
    xml.open(xml_buf);
  }
  catch(const string& e)
  { 
    QDP_error_exit("Error in readMILC: %s",e.c_str());
  }

  END_CODE("readMILC");
}

