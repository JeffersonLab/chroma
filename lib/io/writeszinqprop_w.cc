// $Id: writeszinqprop_w.cc,v 1.1 2003-09-25 22:20:28 edwards Exp $
/*!
 * @file
 * @brief  Write an old SZIN-style (checkerboarded) quark propagator
 */

#include "chromabase.h"
#include "io/writeszinqprop_w.h"

#include "qdp_util.h"   // from QDP++

using namespace QDP;

//! Write a SZIN propagator file. This is a simple memory dump writer.
/*!
 * \ingroup io
 *
 * \param q          propagator ( Read )
 * \param file       path ( Read )
 * \param kappa      kappa value (Read)
 */    

void writeSzinQprop(const LatticePropagator& q, const string& file,
		    const Real& kappa)
{
  BinaryWriter cfg_out(file);

  //
  // Write propagator field
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
