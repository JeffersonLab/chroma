// -*- C++ -*-
// $Id: wallff_w.cc,v 3.0 2006-04-03 04:59:01 edwards Exp $
/*! \file
 *  \brief Structures for wall-sink/source form-factors
 *
 *  Form factors constructed from a quark and a wall sink
 */

#include "chromabase.h"
#include "meas/hadron/wallff_w.h"

namespace Chroma {


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
 * \param t0                 time-slice of the source ( Read )
 */
void wallFormFacSft(multi1d<WallFormFac_momenta_t>& momenta,
		    const LatticeComplex& corr_local_fn,
		    const LatticeComplex& corr_nonlocal_fn,
		    const SftMom& phases,
		    bool compute_nonlocal,
		    int t0)
{
  START_CODE();

  momenta.resize(phases.numMom());  // hold momenta output

  // Length of lattice in decay direction and 3pt correlations fcns
  int length = phases.numSubsets();

  multi2d<DComplex> hsum_local = phases.sft(corr_local_fn);
  multi2d<DComplex> hsum_nonlocal;
  if (compute_nonlocal)
    hsum_nonlocal = phases.sft(corr_nonlocal_fn);
  
  // Loop over insertion momenta
  for(int inser_mom_num=0; inser_mom_num < phases.numMom(); ++inser_mom_num) 
  {
    momenta[inser_mom_num].inser_mom_num = inser_mom_num;
    momenta[inser_mom_num].inser_mom     = phases.numToMom(inser_mom_num);
    
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

    momenta[inser_mom_num].local_current    = local_cur3ptfn;
    momenta[inser_mom_num].nonlocal_current = nonlocal_cur3ptfn;
    
  } // end for(inser_mom_num)

  END_CODE();
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

  write(xml, "gamma_ctr", header.gamma_ctr);
  write(xml, "mu", header.mu);
  write(xml, "gamma_value", header.gamma_value);
  write(xml, "Momenta", header.momenta);

  pop(xml);
}

//! Wallformfac projector writer
void write(XMLWriter& xml, const string& path, const WallFormFac_projector_t& header)
{
  push(xml, path);

  write(xml, "proj_ctr", header.proj_ctr);
  write(xml, "proj_name", header.proj_name);
  write(xml, "Insertion", header.insertion);

  pop(xml);
}

//! Wallformfac lorentz writer
void write(XMLWriter& xml, const string& path, const WallFormFac_lorentz_t& header)
{
  push(xml, path);

  write(xml, "lorentz_ctr", header.lorentz_ctr);
  write(xml, "snk_gamma", header.snk_gamma);
  write(xml, "src_gamma", header.src_gamma);
  write(xml, "Projector", header.projector);

  pop(xml);
}

//! Wallformfac formfac writer
void write(XMLWriter& xml, const string& path, const WallFormFac_formfac_t& header)
{
  push(xml, path);

  write(xml, "formfac_ctr", header.formfac_ctr);
  write(xml, "formfac_name", header.formfac_name);
  write(xml, "Lorentz", header.lorentz);

  pop(xml);
}

//! Wallformfac quark writer
void write(XMLWriter& xml, const string& path, const WallFormFac_quark_t& header)
{
  push(xml, path);

  write(xml, "quark_ctr", header.quark_ctr);
  write(xml, "quark_name", header.quark_name);
  write(xml, "FormFac", header.formfac);

  pop(xml);
}

//! WallFormFac writer
void write(XMLWriter& xml, const string& path, const WallFormFac_formfacs_t& header)
{
  push(xml, path);

  write(xml, "subroutine", header.subroutine);
  write(xml, "Quark", header.quark);

  pop(xml);
}


//! Wallformfac momenta writer
void write(BinaryWriter& bin, const WallFormFac_momenta_t& header)
{
  int magic = 20301;
  write(bin, magic);
  write(bin, header.inser_mom_num);
  write(bin, header.inser_mom);
  write(bin, header.local_current);
  write(bin, header.nonlocal_current);
}

//! Wallformfac insertion writer
void write(BinaryWriter& bin, const WallFormFac_insertion_t& header)
{
  write(bin, header.gamma_ctr);
  write(bin, header.mu);
  write(bin, header.gamma_value);
  write(bin, header.momenta);
}

//! Wallformfac projector writer
void write(BinaryWriter& bin, const WallFormFac_projector_t& header)
{
  write(bin, header.proj_ctr);
  write(bin, header.proj_name);
  write(bin, header.insertion);
}

//! Wallformfac lorentz writer
void write(BinaryWriter& bin, const WallFormFac_lorentz_t& header)
{
  write(bin, header.lorentz_ctr);
  write(bin, header.snk_gamma);
  write(bin, header.src_gamma);
  write(bin, header.projector);
}

//! Wallformfac formfac writer
void write(BinaryWriter& bin, const WallFormFac_formfac_t& header)
{
  write(bin, header.formfac_ctr);
  write(bin, header.formfac_name);
  write(bin, header.lorentz);
}

//! Wallformfac quark writer
void write(BinaryWriter& bin, const WallFormFac_quark_t& header)
{
  write(bin, header.quark_ctr);
  write(bin, header.quark_name);
  write(bin, header.formfac);
}

//! WallFormFac writer
void write(BinaryWriter& bin, const WallFormFac_formfacs_t& header)
{
  write(bin, header.subroutine);
  write(bin, header.quark);
}

}  // end namespace Chroma
