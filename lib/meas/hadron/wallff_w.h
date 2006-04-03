// -*- C++ -*-
// $Id: wallff_w.h,v 3.0 2006-04-03 04:59:01 edwards Exp $
/*! \file
 *  \brief Structures for wall-sink/source form-factors
 *
 *  Form factors constructed from a quark and a wall sink
 */

#ifndef __wallff_h__
#define __wallff_h__

#include "util/ft/sftmom.h"

namespace Chroma {

//! Structures to hold form-factors
struct WallFormFac_momenta_t
{
  int              inser_mom_num;
  multi1d<int>     inser_mom;
  multi1d<Complex> local_current;
  multi1d<Complex> nonlocal_current;
};

struct WallFormFac_insertion_t
{
  int              gamma_ctr;
  int              mu;
  int              gamma_value;
  multi1d<WallFormFac_momenta_t> momenta;
};

struct WallFormFac_projector_t
{
  int              proj_ctr;
  string           proj_name;
  multi1d<WallFormFac_insertion_t>  insertion;
};

struct WallFormFac_lorentz_t
{
  int              lorentz_ctr;
  int              snk_gamma;
  int              src_gamma;
  multi1d<WallFormFac_projector_t>  projector;
};

struct WallFormFac_formfac_t
{
  int              formfac_ctr;
  string           formfac_name;
  multi1d<WallFormFac_lorentz_t>  lorentz;
};

struct WallFormFac_quark_t
{
  int              quark_ctr;
  string           quark_name;
  multi1d<WallFormFac_formfac_t>  formfac;
};

struct WallFormFac_formfacs_t
{
  string           subroutine;
  multi1d<WallFormFac_quark_t>  quark;
};


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
				      const LatticePropagator& anti_prop);


//! Do slow SFT over hadron correlator data
/*!
 * \ingroup hadron
 *
 * \param momenta            momenta structure ( Modify )
 * \param corr_local_fn      contracted local current insertion ( Read )
 * \param corr_nonlocal_fn   contracted nonlocal current insertion ( Read )
 * \param phases             fourier transform phase factors ( Read )
 * \param compute_nonlocal   compute the nonlocal current stuff?? ( Read )
 * \param t0                 time slice of the source ( Read )
 */
void wallFormFacSft(multi1d<WallFormFac_momenta_t>& momenta,
		    const LatticeComplex& corr_local_fn,
		    const LatticeComplex& corr_nonlocal_fn,
		    const SftMom& phases,
		    bool compute_nonlocal,
		    int t0);


// Writers

//! Wallformfac momenta writer
void write(XMLWriter& xml, const string& path, const WallFormFac_momenta_t& header);

//! Wallformfac insertion writer
void write(XMLWriter& xml, const string& path, const WallFormFac_insertion_t& header);

//! Wallformfac projector writer
void write(XMLWriter& xml, const string& path, const WallFormFac_projector_t& header);

//! Wallformfac formfac writer
void write(XMLWriter& xml, const string& path, const WallFormFac_formfac_t& header);

//! Wallformfac lorentz writer
void write(XMLWriter& xml, const string& path, const WallFormFac_lorentz_t& header);

//! Wallformfac quark writer
void write(XMLWriter& xml, const string& path, const WallFormFac_quark_t& header);

//! WallFormFac writer
void write(XMLWriter& xml, const string& path, const WallFormFac_formfacs_t& header);


//! Wallformfac momenta writer
void write(BinaryWriter& bin, const WallFormFac_momenta_t& header);

//! Wallformfac insertion writer
void write(BinaryWriter& bin, const WallFormFac_insertion_t& header);

//! Wallformfac projector writer
void write(BinaryWriter& bin, const WallFormFac_projector_t& header);

//! Wallformfac formfac writer
void write(BinaryWriter& bin, const WallFormFac_formfac_t& header);

//! Wallformfac lorentz writer
void write(BinaryWriter& bin, const WallFormFac_lorentz_t& header);

//! Wallformfac quark writer
void write(BinaryWriter& bin, const WallFormFac_quark_t& header);

//! WallFormFac writer
void write(BinaryWriter& bin, const WallFormFac_formfacs_t& header);

}  // end namespace Chroma

#endif
