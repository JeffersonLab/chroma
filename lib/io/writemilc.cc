// $Id: writemilc.cc,v 1.2 2003-10-14 17:41:23 edwards Exp $

/*! \file
 *  \brief Writer a MILC gauge configuration in the 1997 format
 */

#include "chromabase.h"
#include "io/milc_io.h"
#include "io/writemilc.h"
#include "qdp_util.h"    // from QDP

#include <string>
#include <string.h>
using std::string;

using namespace QDP;

//! Write a MILC configuration file
/*!
 * \ingroup io
 *
 * \param header     structure holding config info ( Read )
 * \param u          gauge configuration ( Read )
 * \param cfg_file   path ( Read )
 */    

void writeMILC(const MILCGauge_t& header, const multi1d<LatticeColorMatrix>& u, 
	       const string& cfg_file)
{
  START_CODE("writeMILC");

  BinaryWriter cfg_out(cfg_file); // for now, cfg_io_location not used

  int magic_number = 20103;
  write(cfg_out, magic_number);

  write(cfg_out, header.nrow, Nd);

  // Time stamp - write exactly 64 bytes padded with nulls
  char date_tmp[65];
  int len = (header.date.size() < 64) ? header.date.size() : 64;
  memset(date_tmp, '\0', 65);
  memcpy(date_tmp, header.date.data(), len);
  cfg_out.writeArray(date_tmp, 64, 1);

  // Site order - only support non-sitelist format
  int order = 0;
  write(cfg_out, order);

  /*
   * Write away...
   */
  
  // MILC format has the directions inside the sites
  multi1d<ColorMatrix> u_old(Nd);

  for(int site=0; site < Layout::vol(); ++site)
  {
    multi1d<int> coord = crtesn(site, Layout::lattSize()); // The coordinate
      
    for(int mu=0; mu < Nd; ++mu)
      u_old[mu] = peekSite(u[mu], coord); // Put it into the correct place

    // Write Nd SU(3) matrices. 
    // NOTE: the su3_matrix layout should be the same as in QDP
    write(cfg_out, u_old, Nd); 
  }

  // Go ahead and write checksums, but will not use for now
  unsigned int sum29=0, sum31=0;  // WARNING: these are BOGUS
  write(cfg_out, sum29);
  write(cfg_out, sum31);

  cfg_out.close();

  END_CODE("writeMILC");
}



//! Write a MILC configuration file
/*!
 * \ingroup io
 *
 * \param xml        xml writer holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Write )
 */    

void writeMILC(XMLBufferWriter& xml, multi1d<LatticeColorMatrix>& u, const string& cfg_file)
{
  START_CODE("writeMILC");

  MILCGauge_t header;
  XMLReader  xml_in(xml);   // use the buffer writer to instantiate a reader
  read(xml_in, "/MILC", header);

  writeMILC(header, u, cfg_file);

  END_CODE("writeMILC");
}

