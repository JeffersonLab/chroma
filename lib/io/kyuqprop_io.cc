// $Id: kyuqprop_io.cc,v 1.1 2004-04-14 20:59:24 edwards Exp $
/*!
 * @file
 * @brief  Read/write a Kentucky quark propagator
 */

#include "chromabase.h"
#include "io/kyuqprop_io.h"
#include "util/ferm/transf.h"

using namespace QDP;

//! Initialize the Dirac to Degrand-Rossi spin transformation matrix
/*!
 * \ingroup io
 *
 * Initialize the similarity transformation matrix from 
 * Euclidean Dirac to Euclidean Degrand-Rossi basis
 *
 * \param U          spin matrix ( Modify )
 *
 * \returns The U in   Gamma_{Degrand-Rossi} = U Gamma_Dirac U^dag
 */    
void initDiracToDRMat(SpinMatrix& U)
{
  /*
   * The magic basis transf is found from
   *
   * NOTE: DR = Degrand-Rossi - the spin basis of QDP
   *
   *  psi_DR = U psi_Dirac
   *  psibar_DR Gamma_DR psi_DR = psibar_Dirac Gamma_Dirac psi_Dirac
   *
   * implies
   *  Gamma_DR = U Gamma_Dirac U^dag
   *
   * and the magic formula is
   *
   *   U = (1/sqrt(2)) | i*sigma_2    i*sigma_2 |
   *                   | i*sigma_2   -i*sigma_2 |
   *     = (1/sqrt(2)) |   0   1        0   1   |
   *                   |  -1   0       -1   0   |
   *                   |   0   1        0  -1   |
   *                   |  -1   0        1   0   |
   *
   *   U^dag = -U = U^transpose
   */
  /*
   * NOTE: I do not see some really short combination of 
   * QDP Gamma matrices that can make this beasty, 
   * so I'll just hardwire it...
   */
  U = zero;
  Complex  one = cmplx( Real(1),Real(0));
  Complex mone = cmplx(-Real(1),Real(0));
  pokeSpin(U,  one, 0, 1);
  pokeSpin(U,  one, 0, 3);
  pokeSpin(U, mone, 1, 0);
  pokeSpin(U, mone, 1, 2);
  pokeSpin(U,  one, 2, 1);
  pokeSpin(U, mone, 2, 3);
  pokeSpin(U, mone, 3, 0);
  pokeSpin(U,  one, 3, 2);

}


//! Read a Kentucky quark propagator
/*!
 * \ingroup io
 *
 * \param q          propagator ( Modify )
 * \param file       path ( Read )
 */    

void readKYUQprop(LatticePropagator& q, const string& file)
{
  BinaryReader cfg_in(file);

  /* KY Indices: 
     x,y,z,t,snk_col,snk_spin,ri,src_col,src_spin 
     x is fastest (Fortran Order)
  */
  LatticePropagator q_old;

  for(int src_spin=0; src_spin < 4; ++src_spin)
    for(int src_color=0; src_color < 3; ++src_color)
    {
//      LatticeReal64  tmp;  // KYU always uses 64 bits
      LatticeDouble  tmp;  // KYU always uses 64 bits

      multi2d<LatticeReal> re(3,4);
      for(int snk_spin=0; snk_spin < 4; ++snk_spin)
	for(int snk_color=0; snk_color < 3; ++snk_color)
	{
	  read(cfg_in, tmp);
	  re(snk_color,snk_spin) = tmp;
	}

      multi2d<LatticeReal> im(3,4);
      for(int snk_spin=0; snk_spin < 4; ++snk_spin)
	for(int snk_color=0; snk_color < 3; ++snk_color)
	{
	  read(cfg_in, tmp);
	  im(snk_color,snk_spin) = tmp;
	}

      // Stuff into a fermion
      LatticeFermion f;
      for(int snk_spin=0; snk_spin < 4; ++snk_spin)
      {
	LatticeColorVector  cv;
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

  cfg_in.close();

  // Now that we have read the prop, need to change the spin basis
  SpinMatrix U;
  initDiracToDRMat(U);

  // And finally...
  q = U * q_old * adj(U);   // note, adj(U) = -U
  
  
}
