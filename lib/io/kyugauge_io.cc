// $Id: kyugauge_io.cc,v 3.1 2007-06-10 14:40:23 edwards Exp $

/*! \file
 *  \brief Read/write a Kentucky gauge configuration
 */

#include "chromabase.h"
#include "io/kyugauge_io.h"

namespace Chroma {


//! Read a Kentucky gauge configuration
/*!
 * \ingroup io
 *
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readKYU(multi1d<LatticeColorMatrix>& u, const string& cfg_file)
{
  START_CODE();

  if (Nc != 3)
  {
    QDPIO::cerr << "readKYU - only supports Nc=3" << endl;
    QDP_abort(1);
  }

  BinaryFileReader cfg_in(cfg_file);

  /* According to Shao Jing the UK config format is:

     u( nxyzt, nri, nc, nc, nd )

     in the Fortran way -- so that means nxyzt changes fastest.

     and nxyzt goes as x + (y-1)*nx + ...
     so that x changes fastest than y than z than t.

     The words are d.p. -- 8 bytes -- or REAL64
  */
  LatticeRealD re, im;
  
  for(int mu=0; mu < Nd; ++mu)
    for(int col=0; col < 3; ++col)
      for(int row=0; row < 3; ++row)
      {
	read(cfg_in, re);
	read(cfg_in, im);

	pokeColor(u[mu], 
		  cmplx(LatticeReal(re),LatticeReal(im)), 
		  row, col);
      }

  cfg_in.close();

  END_CODE();
}


//! Write a Kentucky gauge configuration
/*!
 * \ingroup io
 *
 * \param u          gauge configuration ( Read )
 * \param cfg_file   path ( Read )
 */    

void writeKYU(const multi1d<LatticeColorMatrix>& u, const string& cfg_file)
{
  START_CODE();

  if (Nc != 3)
  {
    QDPIO::cerr << "writeKYU - only supports Nc=3" << endl;
    QDP_abort(1);
  }

  BinaryFileWriter cfg_out(cfg_file);

  /* According to Shao Jing the UK config format is:

     u( nxyzt, nri, nc, nc, nd )

     in the Fortran way -- so that means nxyzt changes fastest.

     and nxyzt goes as x + (y-1)*nx + ...
     so that x changes fastest than y than z than t.

     The words are d.p. -- 8 bytes -- or REAL64
  */
  LatticeComplexD  lc;
  LatticeRealD re, im;
  
  for(int mu=0; mu < Nd; ++mu)
    for(int col=0; col < 3; ++col)
      for(int row=0; row < 3; ++row)
      {
	lc = peekColor(u[mu], row, col);
	re = real(lc);
	im = imag(lc);

	write(cfg_out, re);
	write(cfg_out, im);
      }

  cfg_out.close();

  END_CODE();
}


}  // end namespace Chroma
