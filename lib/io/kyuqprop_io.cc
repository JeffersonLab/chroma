// $Id: kyuqprop_io.cc,v 1.6 2004-05-23 21:44:40 edwards Exp $
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

  if (Ns != 4)
  {
    QDPIO::cerr << "readKYUQprop - only supports Ns=4" << endl;
    QDP_abort(1);
  }

  BinaryReader bin(file);

  /* KY Indices: 
     x,y,z,t,snk_col,snk_spin,ri,src_col,src_spin 
     x is fastest (Fortran Order)
  */
  // KYU always uses 64 bits
  LatticePropagatorD q_old;
  multi2d<LatticeFermionD> source(Nc,Ns);
  multi2d<LatticeRealD> re(Nc,Ns);
  multi2d<LatticeRealD> im(Nc,Ns);
  LatticeFermionD f;
  LatticeColorVectorD  cv;

  for(int src_spin=0; src_spin < Ns; ++src_spin)
    for(int src_color=0; src_color < Nc; ++src_color)
    {
      for(int snk_spin=0; snk_spin < Ns; ++snk_spin)
	for(int snk_color=0; snk_color < Nc; ++snk_color)
	{
	  read(bin, re(snk_color,snk_spin));
	}

      for(int snk_spin=0; snk_spin < Ns; ++snk_spin)
	for(int snk_color=0; snk_color < Nc; ++snk_color)
	{
	  read(bin, im(snk_color,snk_spin));
	}

      // Stuff into a fermion
      for(int snk_spin=0; snk_spin < Ns; ++snk_spin)
      {
	for(int snk_color=0; snk_color < Nc; ++snk_color)
	{
	  pokeColor(cv, 
		    cmplx(re(snk_color,snk_spin), im(snk_color,snk_spin)),
		    snk_color);
	}

	pokeSpin(f, cv, snk_spin);
      }

      // Hold temporarily in a multi2d - will need to rearrange src_spin later
      source(src_color,src_spin) = f;
    }

  bin.close();

  // Now rotate the src_spin indices back to original Dirac basis
  // The source spin had to be in a chiral-like basis to use the 
  // source chirality trick of overlap.
  // Finally, stuff the rotated result into the propagator
  for(int src_color=0; src_color < Nc; ++src_color)
  {
    FermToProp(LatticeFermionD(source(src_color,0)+source(src_color,2)),
	       q_old, src_color, 0);
    FermToProp(LatticeFermionD(source(src_color,1)+source(src_color,3)),
	       q_old, src_color, 1);
    FermToProp(LatticeFermionD(source(src_color,2)-source(src_color,0)),
	       q_old, src_color, 2);
    FermToProp(LatticeFermionD(source(src_color,3)-source(src_color,1)),
	       q_old, src_color, 3);
  }

  // Rescale
  q_old *= RealD(1.0) / sqrt(RealD(2));

  // Now that we have read the prop, need to change the spin basis
  SpinMatrixD U = DiracToDRMat();

  // And finally...
  q = U * q_old * adj(U);   // note, adj(U) = -U
  
  END_CODE("readKYUQprop");
}
