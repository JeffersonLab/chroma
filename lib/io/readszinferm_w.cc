// $Id: readszinferm_w.cc,v 1.1 2003-05-09 20:27:11 edwards Exp $: readszinqprop_w.cc,v 1.6 2003/04/30 21:19:33 edwards Exp $
/*!
 * @file
 * @brief  Read an old SZIN-style (checkerboarded) lattice Dirac fermion
 */

#include "chromabase.h"
#include "io/readszinferm_w.h"

#include "qdp_util.h"   // from QDP++

using namespace QDP;

//! Read a SZIN fermion. This is a simple memory dump reader.
/*!
 * \ingroup io
 *
 * \param q          lattice fermion ( Modify )
 * \param file       path ( Read )
 */    

void readSzinFerm(LatticeFermion& q, const string& file)
{
  BinaryReader cfg_in(file);

  //
  // Read propagator field
  //
  multi1d<int> lattsize_cb = Layout::lattSize();
  Fermion   q_old;

  lattsize_cb[0] /= 2;  // checkerboard in the x-direction in szin

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

      pokeSite(q, q_old, coord); // Put it into the correct place
    }
  }

  cfg_in.close();
}

QDP_END_NAMESPACE();
