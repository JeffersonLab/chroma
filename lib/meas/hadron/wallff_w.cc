// -*- C++ -*-
// $Id: wallff_w.cc,v 1.1 2004-06-02 02:00:03 edwards Exp $
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

