// $Id: unprec_dwf_linop_array_w.cc,v 3.0 2006-04-03 04:58:52 edwards Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_dwf_linop_array_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  //! Creation routine
  /*! \ingroup fermact
   *
   * \param u_            gauge field   (Read)
   * \param WilsonMass_   DWF height    (Read)
   * \param m_q_          quark mass    (Read)
   */
  UnprecDWLinOpArray::UnprecDWLinOpArray(
    Handle< FermState<T,P,Q> > fs,
    const Real& WilsonMass_, const Real& m_q_, int N5_,
    const AnisoParam_t& aniso)
  {
    WilsonMass = WilsonMass_;
    m_q = m_q_;
    a5  = 1.0;
    N5  = N5_;

    fbc = fs->getFermBC();
    D.create(fs, aniso);   // construct using possibly aniso glue

    Real ff = where(aniso.anisoP, aniso.nu / aniso.xi_0, Real(1));
    fact1 =  1 + a5*(1 + (Nd-1)*ff - WilsonMass);
    fact2 = -0.5*a5;
  }



  //! Apply unpreconditioned domain-wall fermion linear operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice
   *
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void UnprecDWLinOpArray::operator() (multi1d<LatticeFermion>& chi, 
				       const multi1d<LatticeFermion>& psi, 
				       enum PlusMinus isign) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    //
    //  Chi   =  D' Psi
    //
    LatticeFermion  tmp;   moveToFastMemoryHint(tmp);

    switch (isign)
    {
    case PLUS:
      for(int n=0; n < N5; ++n)
      {
	D(tmp, psi[n], isign);

	if (n == 0)
	  chi[n] = fact2*tmp + fact1*psi[n] 
	    + m_q*chiralProjectPlus(psi[N5-1]) - chiralProjectMinus(psi[1]);
	else if (n == N5-1)
	  chi[n] = fact2*tmp + fact1*psi[n] 
	    - chiralProjectPlus(psi[N5-2]) + m_q*chiralProjectMinus(psi[0]);
	else
	  chi[n] = fact2*tmp + fact1*psi[n] 
	    - chiralProjectPlus(psi[n-1]) - chiralProjectMinus(psi[n+1]);
      }          
      break;

    case MINUS:
      for(int n=0; n < N5; ++n)
      {
	D(tmp, psi[n], isign);

	if (n == 0)
	  chi[n] = fact2*tmp + fact1*psi[n] 
	    + m_q*chiralProjectMinus(psi[N5-1]) - chiralProjectPlus(psi[1]);
	else if (n == N5-1)
	  chi[n] = fact2*tmp + fact1*psi[n] 
	    - chiralProjectMinus(psi[N5-2]) + m_q*chiralProjectPlus(psi[0]);
	else
	  chi[n] = fact2*tmp + fact1*psi[n] 
	    - chiralProjectMinus(psi[n-1]) - chiralProjectPlus(psi[n+1]);
      }          
      break;
    }

    getFermBC().modifyF(chi);

    END_CODE();
  }

  //! Derivative
  void 
  UnprecDWLinOpArray::deriv(multi1d<LatticeColorMatrix>& ds_u, 
			    const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
			    enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);
    ds_u = zero;

    Real fact2 = -Real(0.5)*a5;

    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    for(int s(0);s<N5;s++)
    {
      D.deriv(ds_tmp,chi[s],psi[s],isign);
      for(int mu(0);mu<Nd;mu++)
	ds_u[mu] += fact2*ds_tmp[mu];
    }

    END_CODE();
  }

}; // End Namespace Chroma

