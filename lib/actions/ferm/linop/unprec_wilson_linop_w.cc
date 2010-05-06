// $Id: unprec_wilson_linop_w.cc,v 3.2 2006-08-26 05:50:06 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "handle.h"
#include "state.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 

  //! Creation routine
  /*!
   * \param u_ 	  gauge field     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   */
  void UnprecWilsonLinOp::create(Handle< FermState<T,P,Q> > fs,
				 const Real& Mass_)
  {
    AnisoParam_t anisoParam;
    create(fs, Mass_, anisoParam);
  }


  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	  gauge field     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   * \param aniso   anisotropy struct   	       (Read)
   */
  void UnprecWilsonLinOp::create(Handle< FermState<T,P,Q> > fs,
				 const Real& Mass_,
				 const AnisoParam_t& anisoParam)
  {
    START_CODE();

    D.create(fs,anisoParam);

    Mass = Mass_;
    Real ff = where(anisoParam.anisoP, anisoParam.nu / anisoParam.xi_0, Real(1));
    fact = 1 + (Nd-1)*ff + Mass;
    
    END_CODE();
  }


  //! Apply unpreconditioned Wilson fermion linear operator
  /*!
   * \ingroup linop
   *
   * The operator acts on the entire lattice
   *
   * \param chi 	  Pseudofermion field     	       (Read)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void UnprecWilsonLinOp::operator() (LatticeFermion& chi, const LatticeFermion& psi, 
				      enum PlusMinus isign) const
  {
    START_CODE();

    //
    //  Chi   =  (Nd+Mass)*Psi  -  (1/2) * D' Psi
    //
    LatticeFermion tmp;   moveToFastMemoryHint(tmp);
    Real mhalf = -0.5;

    // D is a Dslash - must apply to both CB-s
    D(tmp, psi, isign);

    chi = fact*psi + mhalf*tmp;

    getFermBC().modifyF(chi);
  
    END_CODE();
  }

    //! Apply operator with towers
  void UnprecWilsonLinOp::operator()(Tower<T>& chi, const Tower<T>& psi,
				     const P& p,
				     enum PlusMinus isign)
  {
 
    if( chi.size() != psi.size() ) { 
      QDPIO::cout << "Error towers are incompatible" << endl;
      QDP_abort(1);
    }
    int N = chi.size();
    Real mhalf = -0.5;

    D(chi, psi, p, isign);   // All the fun stuff is in here.
    chi *= mhalf;

    // This is an optimized tower addition. Really it corresponds to
    // [ 0, 0,..., fact ] * [ psi_N, psi_N-1, ..., psi_0 ]
    // = [fact*psi_N, fact*psi_N-1,..., fact*psi_N ]
    //
    for(int i=0; i < N; i++) {
      chi[i] += fact*psi[i];
    }
    
    getFermBC().modifyF(chi);
  
    END_CODE();
  }

  void UnprecWilsonLinOp::deriv(TowerArray<PQTraits<Q>::Base_t>& ds_u,
	       const Tower<T>& chi,
	       const Tower<T>& psi,
	       const P& p,
	       enum PlusMinus isign)
    {
      
      // This does both parities
      D.deriv(ds_u, chi, psi, isign);

      // Factor from the -1/2 in front of the dslash
      for(int mu = 0; mu < Nd; ++mu)
	ds_u[mu] *= Real(-0.5);


    }


  //! Derivative of unpreconditioned Wilson dM/dU
  /*!
   * \param chi     left vector on cb                           (Read)
   * \param psi     right vector on 1-cb                        (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb	    Checkerboard of chi vector                  (Read)
   *
   * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
   */
  void 
  UnprecWilsonLinOp::deriv(multi1d<LatticeColorMatrix>& ds_u,
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
  {
    START_CODE();
    // This does both parities
    D.deriv(ds_u, chi, psi, isign);

    // Factor from the -1/2 in front of the dslash
    for(int mu = 0; mu < Nd; ++mu)
      ds_u[mu] *= Real(-0.5);

    END_CODE();
  }


  //! Return flops performed by the operator()
  unsigned long UnprecWilsonLinOp::nFlops() const
  {
    unsigned long site_flops = D.nFlops()+4*Nc*Ns;
    return site_flops*Layout::sitesOnNode();
  }

} // End Namespace Chroma
