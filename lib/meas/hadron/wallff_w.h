// -*- C++ -*-
// $Id: wallff_w.h,v 1.1 2004-06-02 02:00:03 edwards Exp $
/*! \file
 *  \brief Structures for wall-sink/source form-factors
 *
 *  Form factors constructed from a quark and a wall sink
 */

#ifndef __wallff_h__
#define __wallff_h__


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
  int              gamma_value;
  multi1d<WallFormFac_momenta_t> momenta;
};

struct WallFormFac_insertions_t
{
  int              seq_src;
  multi1d<WallFormFac_insertion_t>  insertions;
};

struct WallFormFac_formfacs_t
{
  multi1d<WallFormFac_insertions_t>  formFacs;
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


// Writers

//! Wallformfac momenta writer
void write(XMLWriter& xml, const string& path, const WallFormFac_momenta_t& header);

//! Wallformfac insertion writer
void write(XMLWriter& xml, const string& path, const WallFormFac_insertion_t& header);

//! Wallformfac insertions writer
void write(XMLWriter& xml, const string& path, const WallFormFac_insertions_t& header);

//! WallFormFac writer
void write(XMLWriter& xml, const string& path, const WallFormFac_formfacs_t& header);


#endif
