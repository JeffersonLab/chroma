// $Id: readszinqprop_w.cc,v 1.10 2003-10-10 03:46:46 edwards Exp $
/*!
 * @file
 * @brief  Read an old SZIN-style (checkerboarded) quark propagator
 */

#include "chromabase.h"
#include "io/readszinqprop_w.h"

#include "qdp_util.h"   // from QDP++

using namespace QDP;

//! Read a SZIN propagator file. This is a simple memory dump reader.
/*!
 * \ingroup io
 *
 * \param xml        xml reader holding prop info ( Modify )
 * \param q          propagator ( Modify )
 * \param file       path ( Read )
 */    

void readSzinQprop(XMLReader& xml, LatticePropagator& q, const string& file)
{
  BinaryReader cfg_in(file);

  //
  // Read propagator field
  //
  multi1d<int> lattsize_cb = Layout::lattSize();
  Propagator   q_tmp, q_old;
  Real Kappa;

  lattsize_cb[0] /= 2;  // checkerboard in the x-direction in szin

  // Read Kappa
  read(cfg_in, Kappa);

  // Read prop
  for(int cb=0; cb < 2; ++cb)
  {
    for(int sitecb=0; sitecb < Layout::vol()/2; ++sitecb)
    {
      multi1d<int> coord = crtesn(sitecb, lattsize_cb);

      // construct the checkerboard offset
      int sum = 0;
      for(int m=1; m<Nd; m++)
	sum += coord[m];

      // The true lattice x-coord
      coord[0] = 2*coord[0] + ((sum + cb) & 1);

      read(cfg_in,q_old); 	// Read in a site propagator

      q_tmp = transpose(q_old); // Take the transpose in both color and spin space
      pokeSite(q, q_tmp, coord); // Put it into the correct place
    }
  }

  cfg_in.close();


  // Now, set up the XML header. Do this by first making a buffer
  // writer that is then used to make the reader
  XMLBufferWriter  xml_buf;

  push(xml_buf, "szin_prop");
  write(xml_buf,"Kappa",Kappa);
  pop(xml_buf);

  try 
  {
    xml.open(xml_buf);
  }
  catch(const string& e)
  { 
    QDP_error_exit("Error in readszinqprop: %s",e.c_str());
  }

}

//! Write a SZIN propagator file. This is a simple memory dump writer.
/*!
 * \ingroup io
 *
 * \param q          propagator ( Read )
 * \param file       path ( Read )
 * \param kappa      kappa value (Read)
 */    

void writeSzinQprop(const LatticePropagator& q, const string& file,
		    const Real kappa)
{
  BinaryWriter cfg_out(file);

  //
  // Read propagator field
  //
  multi1d<int> lattsize_cb = Layout::lattSize();
  Propagator   q_tmp, q_old;

  lattsize_cb[0] /= 2;  // checkerboard in the x-direction in szin

  // Write Kappa
  write(cfg_out, kappa);

  // Write prop
  for(int cb=0; cb < 2; ++cb)
  {
    for(int sitecb=0; sitecb < Layout::vol()/2; ++sitecb)
    {
      multi1d<int> coord = crtesn(sitecb, lattsize_cb);

      // construct the checkerboard offset
      int sum = 0;
      for(int m=1; m<Nd; m++)
	sum += coord[m];

      // The true lattice x-coord
      coord[0] = 2*coord[0] + ((sum + cb) & 1);

      q_old = peekSite(q, coord); // Get the correct site
      
      q_tmp = transpose(q_old);	// Take the transpose
      write(cfg_out, q_tmp);

    }
  }

  cfg_out.close();
}
