// $Id: readszinqprop_w.cc,v 3.1 2007-06-10 14:40:23 edwards Exp $
/*!
 * @file
 * @brief  Read an old SZIN-style (checkerboarded) quark propagator
 */

#include <iostream>
#include <string>
#include <sstream>

#include "chromabase.h"
#include "io/readszinqprop_w.h"

#include "qdp_util.h"   // from QDP++

namespace Chroma {

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
  BinaryFileReader cfg_in(file);

  //
  // Read propagator field
  //
  multi1d<int> lattsize_cb = Layout::lattSize();
  Real Kappa;

  lattsize_cb[0] /= 2;  // checkerboard in the x-direction in szin

  // Read Kappa
  read(cfg_in, Kappa);

  // Read prop
  LatticePropagator  q_old;

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

      read(cfg_in, q_old, coord); 	// Read in a site propagator
    }
  }

  q = transpose(q_old);  // take the transpose

  cfg_in.close();


  // Now, set up the XML header. Do this by first making a buffer
  // writer that is then used to make the reader
  XMLBufferWriter  xml_buf;

  push(xml_buf, "szin_prop");
  write(xml_buf,"Kappa",Kappa);
  pop(xml_buf);

  try 
  {
    // Temporary XLC failure workaround
#if 0
    xml.open(xml_buf);
#else
    const string bufcontent=xml_buf.str() + "\n";
    istringstream is(bufcontent);
    xml.open(is);
#endif

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
  BinaryFileWriter cfg_out(file);

  //
  // Read propagator field
  //
  multi1d<int> lattsize_cb = Layout::lattSize();

  lattsize_cb[0] /= 2;  // checkerboard in the x-direction in szin

  // Write Kappa
  write(cfg_out, kappa);

  // Write prop
  LatticePropagator  q_old = transpose(q);   // take the transpose

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

      write(cfg_out, q_old, coord);   // write out a single site
    }
  }

  cfg_out.close();
}

}  // end namespace Chroma
