// $Id: kyugauge_io.cc,v 1.3 2004-04-15 03:33:54 edwards Exp $

/*! \file
 *  \brief Read a Kentucky gauge configuration
 */

#include "chromabase.h"
#include "io/kyugauge_io.h"

using namespace QDP;

//! Read a Kentucky gauge configuration
/*!
 * \ingroup io
 *
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readKYU(multi1d<LatticeColorMatrix>& u, const string& cfg_file)
{
  START_CODE("readKYU");

  if (Nc != 3)
  {
    QDPIO::cerr << "readKYU - only supports Nc=3" << endl;
    QDP_abort(1);
  }

  BinaryReader cfg_in(cfg_file); // for now, cfg_io_location not used

  /* According to Shao Jing the UK config format is:

     u( nxyzt, nri, nc, nc, nd )

     in the Fortran way -- so that means nxyzt changes fastest.

     and nxyzt goes as x + (y-1)*nx + ...
     so that x changes fastest than y than z than t.

     The words are d.p. -- 8 bytes -- or REAL64 
  */
//  LatticeReal64 re, im;
  LatticeDouble re, im;
  
  for(int mu=0; mu < Nd; ++mu)
    for(int col=0; col < 3; ++col)
      for(int row=0; row < 3; ++row)
      {
	read(cfg_in, re);
	read(cfg_in, im);

	pokeColor(u[mu], 
		  cmplx(LatticeReal(re),LatticeReal(im)), 
		  col, row);   // transpose
      }

  cfg_in.close();

  END_CODE("readKYU");
}

