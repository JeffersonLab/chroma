// -*- C++ -*-
// $Id: wallff_w.cc,v 1.2 2004-06-04 04:05:58 edwards Exp $
/*! \file
 *  \brief Structures for wall-sink/source form-factors
 *
 *  Form factors constructed from a quark and a wall sink
 */

#include "chromabase.h"
#include "meas/hadron/wallff_w.h"

using namespace QDP;


//! Compute nonlocal current propagator
/*!
 * \ingroup hadron
 *
 * The form of J_mu = (1/2)*[psibar(x+mu)*U^dag_mu*(1+gamma_mu)*psi(x) -
 *                           psibar(x)*U_mu*(1-gamma_mu)*psi(x+mu)]
 *
 * \param u                  gauge fields ( Read )
 * \param mu                 direction ( Read )
 * \param forw_prop          forward propagator ( Read )
 * \param anti_prop          anti-quark version of forward propagator ( Read )
 *
 * \return nonlocal current propagator
 */
LatticePropagator nonlocalCurrentProp(const multi1d<LatticeColorMatrix>& u, 
				      int mu, 
				      const LatticePropagator& forw_prop,
				      const LatticePropagator& anti_prop)
{
  int gamma_value = 1 << mu;

  LatticePropagator S = 0.5*(shift(anti_prop, FORWARD, mu) * adj(u[mu])
    * (forw_prop + Gamma(gamma_value)*forw_prop)
    - anti_prop * u[mu] * shift(forw_prop - Gamma(gamma_value)*forw_prop, FORWARD, mu));

  return S;
}


//! Do slow SFT over hadron correlator data
/*!
 * \ingroup hadron
 *
 * \param momenta            momenta structure ( Modify )
 * \param corr_local_fn      contracted local current insertion ( Read )
 * \param corr_nonlocal_fn   contracted nonlocal current insertion ( Read )
 * \param phases             fourier transform phase factors ( Read )
 * \param compute_nonlocal   compute the nonlocal current stuff?? ( Read )
 * \param t0                 time coordinates of the source ( Read )
 * \param t_sink             time coordinates of the sink ( Read )
 */
void wallFormFacSft(multi1d<WallFormFac_momenta_t>& momenta,
		    const LatticeComplex& corr_local_fn,
		    const LatticeComplex& corr_nonlocal_fn,
		    const SftMom& phases,
		    bool compute_nonlocal,
		    int t0, int t_sink)
{
  START_CODE("wallFormFacSft");

  if ( momenta.size() != phases.numMom() )
  {
    QDPIO::cerr << "wallFormFacSft: momenta array size incorrect" << endl;
    QDP_abort(1);
  }

  // Length of lattice in j_decay direction and 3pt correlations fcns
  int length = phases.numSubsets();

  multi2d<DComplex> hsum_local = phases.sft(corr_local_fn);
  multi2d<DComplex> hsum_nonlocal;
  if (compute_nonlocal)
    hsum_nonlocal = phases.sft(corr_nonlocal_fn);
  
  // Loop over insertion momenta
  for(int inser_mom_num=0; inser_mom_num < phases.numMom(); ++inser_mom_num) 
  {
    momenta[insert_mom_num].inser_mom_num = inser_mom_num;
    momenta[insert_mom_num].inser_mom     = phases.numToMom(inser_mom_num);
    
    multi1d<Complex> local_cur3ptfn(length); // always compute
    multi1d<Complex> nonlocal_cur3ptfn;
    if (compute_nonlocal)
      nonlocal_cur3ptfn.resize(length);      // possibly compute
	    
    for (int t=0; t < length; ++t) 
    {
      int t_eff = (t - t0 + length) % length;
      
      local_cur3ptfn[t_eff] = Complex(hsum_local[inser_mom_num][t]);
      if (compute_nonlocal)
	nonlocal_cur3ptfn[t_eff] = Complex(hsum_nonlocal[inser_mom_num][t]);
    } // end for(t)

    momenta[insert_mom_num].local_current    = local_cur3ptfn;
    momenta[insert_mom_num].nonlocal_current = nonlocal_cur3ptfn;
    
  } // end for(inser_mom_num)

  END_CODE("wallFormFacSft");
}




// Writers
//! Wallformfac momenta writer
void write(XMLWriter& xml, const string& path, const WallFormFac_momenta_t& header)
{
  push(xml, path);

  write(xml, "inser_mom_num", header.inser_mom_num);
  write(xml, "inser_mom", header.inser_mom);
  write(xml, "local_cur3ptfn", header.local_current);

  if (header.nonlocal_current.size() > 0)
    write(xml, "nonlocal_cur3ptfn", header.nonlocal_current);

  pop(xml);
}

//! Wallformfac insertion writer
void write(XMLWriter& xml, const string& path, const WallFormFac_insertion_t& header)
{
  push(xml, path);

  write(xml, "gamma_value", header.gamma_value);
  write(xml, "Momenta", header.momenta);

  pop(xml);
}

//! Wallformfac insertions writer
void write(XMLWriter& xml, const string& path, const WallFormFac_insertions_t& header)
{
  push(xml, path);

  write(xml, "seq_src", header.seq_src);
  write(xml, "Insertions", header.insertions);

  pop(xml);
}

//! WallFormFac writer
void write(XMLWriter& xml, const string& path, const WallFormFac_formfacs_t& header)
{
  write(xml, path, header.formFacs);
}

