// $Id: kyuqprop_io.cc,v 1.4 2004-05-14 00:27:35 edwards Exp $
/*!
 * @file
 * @brief  Read/write a Kentucky quark propagator
 */

#include "chromabase.h"
#include "io/kyuqprop_io.h"
#include "util/ferm/transf.h"
#include "util/ferm/diractodr.h"

using namespace QDP;

//! Read a Kentucky quark propagator
/*!
 * \ingroup io
 *
 * \param q          propagator ( Modify )
 * \param file       path ( Read )
 */    

void readKYUQprop(LatticePropagator& q, const string& file)
{
  START_CODE("readKYUQprop");

  if (Nc != 3)
  {
    QDPIO::cerr << "readKYUQprop - only supports Nc=3" << endl;
    QDP_abort(1);
  }

  BinaryReader bin(file);

  /* KY Indices: 
     x,y,z,t,snk_col,snk_spin,ri,src_col,src_spin 
     x is fastest (Fortran Order)
  */
  LatticePropagator q_old;
  multi2d<LatticeReal> re(3,4);
  multi2d<LatticeReal> im(3,4);
  LatticeFermion f;
  LatticeColorVector  cv;

//  LatticeReal64  tmp;  // KYU always uses 64 bits
  LatticeDouble  tmp;  // KYU always uses 64 bits

  for(int src_spin=0; src_spin < 4; ++src_spin)
    for(int src_color=0; src_color < 3; ++src_color)
    {
      for(int snk_spin=0; snk_spin < 4; ++snk_spin)
	for(int snk_color=0; snk_color < 3; ++snk_color)
	{
	  read(bin, tmp);
	  re(snk_color,snk_spin) = tmp;
	}

      for(int snk_spin=0; snk_spin < 4; ++snk_spin)
	for(int snk_color=0; snk_color < 3; ++snk_color)
	{
	  read(bin, tmp);
	  im(snk_color,snk_spin) = tmp;
	}

      // Stuff into a fermion
      for(int snk_spin=0; snk_spin < 4; ++snk_spin)
      {
	for(int snk_color=0; snk_color < 3; ++snk_color)
	{
	  pokeColor(cv, 
		    cmplx(re(snk_color,snk_spin), im(snk_color,snk_spin)),
		    snk_color);
	}

	pokeSpin(f, cv, snk_spin);
      }

      // Stuff into the propagator
      FermToProp(f, q_old, src_color, src_spin);
    }

  bin.close();

  // Now that we have read the prop, need to change the spin basis
  SpinMatrix U = DiracToDRMat();

  // And finally...
  q = U * q_old * adj(U);   // note, adj(U) = -U
  
  END_CODE("readKYUQprop");
}
