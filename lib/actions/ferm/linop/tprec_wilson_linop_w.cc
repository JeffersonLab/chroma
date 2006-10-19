// $Id: tprec_wilson_linop_w.cc,v 3.1 2006-10-19 16:01:32 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_wilson_linop_w.h"

using namespace QDP::Hints;
namespace Chroma 
{ 
  //! Creation routine
  /*!
   * \param u_ 	  gauge field     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   */
  void EvenOddPrecWilsonLinOp::create(Handle< FermState<T,P,Q> > fs,
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
  void EvenOddPrecWilsonLinOp::create(Handle< FermState<T,P,Q> > fs,
				      const Real& Mass_,
				      const AnisoParam_t& anisoParam)
  {
    START_CODE();

    D.create(fs,anisoParam);

    Mass = Mass_;
    Real ff = where(anisoParam.anisoP, anisoParam.nu / anisoParam.xi_0, Real(1));
    fact = 1 + (Nd-1)*ff + Mass;
    invfact = 1/fact;
    
    END_CODE();
  }


  //! Apply even-odd linop component
  /*!
   * The operator acts on the entire even sublattice
   *
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void 
  EvenOddPrecWilsonLinOp::evenOddLinOp(LatticeFermion& chi, 
				       const LatticeFermion& psi, 
				       enum PlusMinus isign) const
  {
    START_CODE();

    Real mhalf = -0.5;

    D.apply(chi, psi, isign, 0);
    chi[rb[0]] *= mhalf;
  
    END_CODE();
  }

  //! Apply odd-even linop component
  /*!
   * The operator acts on the entire odd sublattice
   *
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void 
  EvenOddPrecWilsonLinOp::oddEvenLinOp(LatticeFermion& chi, 
				       const LatticeFermion& psi, 
				       enum PlusMinus isign) const
  {
    START_CODE();

    Real mhalf = -0.5;

    D.apply(chi, psi, isign, 1);
    chi[rb[1]] *= mhalf;
  
    END_CODE();
  }



  //! Override inherited one with a few more funkies
  void EvenOddPrecWilsonLinOp::operator()(LatticeFermion & chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
  {
    START_CODE();

    LatticeFermion tmp1, tmp2, tmp3;  // if an array is used here, 

    moveToFastMemoryHint(tmp1);
    moveToFastMemoryHint(tmp2);
    moveToFastMemoryHint(tmp3);
  

    Real mquarterinvfact = -0.25*invfact;

    // tmp1[0] = D_eo psi[1]
    D.apply(tmp1, psi, isign, 0);

    // tmp2[1] = D_oe tmp1[0]
    D.apply(tmp2, tmp1, isign, 1);

    // Now we have tmp2[1] = D_oe D_eo psi[1]

    // now scale tmp2[1] with (-1/4)/fact = (-1/4)*(1/(Nd + m))
    // with a vscale -- using tmp2 on both sides should be OK, but
    // just to be safe use tmp3
    tmp3[rb[1]] = mquarterinvfact*tmp2;

    // now tmp3[1] should be = (-1/4)*(1/(Nd + m) D_oe D_eo psi[1]

    // Now get chi[1] = fact*psi[1] + tmp3[1]
    // with a vaxpy3 
    // chi[1] = (Nd + m) - (1/4)*(1/(Nd + m)) D_oe D_eo psi[1]
    //
    // in this order, this last job could be replaced with a 
    // vaxpy3_norm if we wanted the || M psi ||^2 over the subset.
    chi[rb[1]] = fact*psi + tmp3;

    getFermBC().modifyF(chi, rb[1]);
    
    END_CODE();
  }


  //! Derivative of even-odd linop component
  void 
  EvenOddPrecWilsonLinOp::derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u,
					    const LatticeFermion& chi, 
					    const LatticeFermion& psi, 
					    enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);

    D.deriv(ds_u, chi, psi, isign, 0);
    for(int mu=0; mu < Nd; mu++) {
      ds_u[mu] *=  Real(-0.5);
    }

    END_CODE();
  }


  //! Derivative of odd-even linop component
  void 
  EvenOddPrecWilsonLinOp::derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
					    const LatticeFermion& chi, 
					    const LatticeFermion& psi, 
					    enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);

    D.deriv(ds_u, chi, psi, isign, 1);
    for(int mu=0; mu < Nd; mu++) { 
      ds_u[mu]  *= Real(-0.5);
    }
    END_CODE();
  }

  //! Return flops performed by the operator()
  unsigned long EvenOddPrecWilsonLinOp::nFlops() const
  { 
    unsigned long cbsite_flops = 2*D.nFlops()+6*Nc*Ns;
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

} // End Namespace Chroma
