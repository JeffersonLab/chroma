// $Id: coulgauge.cc,v 3.1 2009-10-09 15:33:43 bjoo Exp $
/*! \file
 *  \brief Coulomb (and Landau) gauge fixing 
 */

#include "chromabase.h"
#include "meas/gfix/coulgauge.h"
#include "meas/gfix/grelax.h"
#include "util/gauge/reunit.h"

namespace Chroma {

/********************** HACK ******************************/
// Primitive way for now to indicate the time direction
static int tDir() {return Nd-1;}
static Real xi_0() {return 1.0;}
/******************** END HACK ***************************/


//! Coulomb (and Landau) gauge fixing
/*!
 * \ingroup gfix
 *
 * Driver for gauge fixing to Coulomb gauge in slices perpendicular
 * to the direction "j_decay".
 * If j_decay >= Nd: fix to Landau gauge.
 * Note: as written this works only for SU(2) and SU(3)!

 * \param u        (gauge fixed) gauge field ( Modify )
 * \param n_gf     number of gauge fixing iterations ( Write )
 * \param j_decay  direction perpendicular to slices to be gauge fixed ( Read )
 * \param GFAccu   desired accuracy for gauge fixing ( Read )
 * \param GFMax    maximal number of gauge fixing iterations ( Read )
 * \param OrDo     use overrelaxation or not ( Read )
 * \param OrPara   overrelaxation parameter ( Read )
 */

void coulGauge(multi1d<LatticeColorMatrix>& u, 
	       int& n_gf, 
	       int j_decay, const Real& GFAccu, int GFMax, 
	       bool OrDo, const Real& OrPara)
{
  LatticeColorMatrix g;

  coulGauge(u, g, n_gf, j_decay, GFAccu, GFMax, OrDo, OrPara);
}



//! Coulomb (and Landau) gauge fixing
/*!
 * \ingroup gfix
 *
 * Driver for gauge fixing to Coulomb gauge in slices perpendicular
 * to the direction "j_decay".
 * If j_decay >= Nd: fix to Landau gauge.
 * Note: as written this works only for SU(2) and SU(3)!

 * \param u        (gauge fixed) gauge field ( Modify )
 * \param g        Gauge transformation matrices (Write)
 * \param n_gf     number of gauge fixing iterations ( Write )
 * \param j_decay  direction perpendicular to slices to be gauge fixed ( Read )
 * \param GFAccu   desired accuracy for gauge fixing ( Read )
 * \param GFMax    maximal number of gauge fixing iterations ( Read )
 * \param OrDo     use overrelaxation or not ( Read )
 * \param OrPara   overrelaxation parameter ( Read )
 */

void coulGauge(multi1d<LatticeColorMatrix>& u, 
	       LatticeColorMatrix& g,
	       int& n_gf, 
	       int j_decay, const Real& GFAccu, int GFMax, 
	       bool OrDo, const Real& OrPara)
{
  Double tgfold;
  Double tgfnew;
  Double tgf_t;
  Double tgf_s;
  Double norm;
  int num_sdir;
  bool tdirp;

  START_CODE();

  Real xi_sq = pow(xi_0(),2);
  if( j_decay >= 0 && j_decay < Nd )
  {
    if( tDir() >= 0 && tDir() != j_decay )
    {
      num_sdir = Nd - 2;
      tdirp = true;
      norm = Double(Layout::vol()*Nc) * (Double(num_sdir)+Double(xi_sq));
    }
    else
    {
      num_sdir = Nd - 1;
      tdirp = false;
      norm = Double(Layout::vol()*Nc*num_sdir);
    }
  }
  else
  {
    if( tDir() >= 0 && tDir() < Nd )
    {
      num_sdir = Nd - 1;
      tdirp = true;
      norm = Double(Layout::vol()*Nc) * (Double(num_sdir)+Double(xi_sq));
    }
    else
    {
      num_sdir = Nd;
      tdirp = false;
      norm = Double(Layout::vol()*Nc*num_sdir);
    }
  }

      
  /* Compute initial gauge fixing term: sum(trace(U_spacelike)); */
  tgf_t = 0;
  tgf_s = 0;
  for(int mu=0; mu<Nd; ++mu)
    if( mu != j_decay )
    {
      Double tgf_tmp = sum(real(trace(u[mu])));

      if( mu != tDir() )
	tgf_s += tgf_tmp;
      else
	tgf_t += tgf_tmp;
    }

  if( tdirp )
  {
    tgfold = (xi_sq*tgf_t+tgf_s)/norm;
    tgf_s = tgf_s/(Double(Layout::vol()*Nc*num_sdir));
    tgf_t = tgf_t/(Double(Layout::vol()*Nc));
  }
  else
  {
    tgf_s = tgf_s/(Double(Layout::vol()*Nc*num_sdir));
    tgfold = tgf_s;
  }
  
  // Gauge transf. matrices always start from identity
  g = 1; 

  /* Gauge fix until converged or too many iterations */
  n_gf = 0;
  bool wrswitch = true;    /* switch for writing of gauge fixing term */
  Double conver = 1;        /* convergence criterion */

  while( toBool(conver > GFAccu)  &&  n_gf < GFMax )
  {
    n_gf = n_gf + 1;
    if( GFMax - n_gf < 11 ) 
      wrswitch = true;
    
    /* Loop over checkerboards for gauge fixing */
    for(int cb=0; cb<2; ++cb)
    {
      if (Nc > 1)
      {
	/* Loop over SU(2) subgroup index */
	for(int su2_index=0; su2_index < Nc*(Nc-1)/2; ++su2_index)
	{
	  /* Now do a gauge fixing relaxation step */
	  grelax(g, u, j_decay, su2_index, cb, OrDo, OrPara);
	}   /* end su2_index loop */
      }
      else
      {
	int su2_index = -1;
	/* Now do a gauge fixing relaxation step */
	grelax(g, u, j_decay, su2_index, cb, OrDo, OrPara);
      }
    }     /* end cb loop */

    /* Reunitarize */
    reunit(g);

    /* Compute new gauge fixing term: sum(trace(U_spacelike)): */
    tgf_t = 0;
    tgf_s = 0;
    for(int mu=0; mu<Nd; ++mu)
      if( mu != j_decay )
      {
	Double tgf_tmp = sum(real(trace(g * u[mu] * shift(adj(g), FORWARD, mu))));

	if( mu != tDir() )
	  tgf_s += tgf_tmp;
	else
	  tgf_t += tgf_tmp;
      }

    if( tdirp )
    {
      tgfnew = (xi_sq*tgf_t+tgf_s)/norm;
      tgf_s = tgf_s/(Double(Layout::vol()*Nc*num_sdir));
      tgf_t = tgf_t/(Double(Layout::vol()*Nc));
    }
    else
    {
      tgf_s = tgf_s/(Double(Layout::vol()*Nc*num_sdir));
      tgfnew = tgf_s;
    }

    if( wrswitch ) 
      QDPIO::cout << "COULGAUGE: iter= " << n_gf 
		  << "  tgfold= " << tgfold 
		  << "  tgfnew= " << tgfnew
		  << "  tgf_s= " << tgf_s 
		  << "  tgf_t= " << tgf_t << endl;

    /* Normalized convergence criterion: */
    conver = fabs((tgfnew - tgfold) / tgfnew);
    tgfold = tgfnew;
  }       /* end while loop */

      
  if( wrswitch )
    QDPIO::cout << "COULGAUGE: end: iter= " << n_gf 
		<< "  tgfold= " << tgfold 
		<< "  tgf_s= " << tgf_s 
		<< "  tgf_t= " << tgf_t << endl;

  // Finally, gauge rotate the original matrices and overwrite them
  for(int mu = 0; mu < Nd; ++mu)
  {
    LatticeColorMatrix u_tmp = g * u[mu];
    u[mu] = u_tmp * shift(adj(g), FORWARD, mu);
  }
    
#if 0
  /*+ debugging */
  XMLBufferWriter xml_out;
  push(xml_out,"Final_trace_max_in_CoulGauge");
  write(xml_out, "j_decay", j_decay);
  write(xml_out, "t_dir", tDir());
  write(xml_out, "n_gf",n_gf);
  write(xml_out, "tgfold", tgfold);
  write(xml_out, "tgf_s", tgf_s);
  write(xml_out, "tgf_t", tgf_t);
  pop(xml_out);
#endif

  END_CODE();
}


}; // Namespace Chroma
