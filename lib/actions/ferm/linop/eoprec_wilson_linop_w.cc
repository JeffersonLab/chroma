// $Id: eoprec_wilson_linop_w.cc,v 3.2 2008-01-21 20:18:50 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson linear operator
 */

#include "actions/ferm/linop/eoprec_wilson_linop_w.h"

using namespace QDP::Hints;
namespace Chroma 
{ 
  //! Creation routine
  /*!
   * \param fs 	    gauge state     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   */
  void EvenOddPrecWilsonLinOp::create(Handle< FermState<T,P,Q> > fs,
				      const Real& Mass_)
  {
    multi1d<Real> cf(Nd);
    cf = 1.0;
    create(fs, Mass_, cf);
  }


  //! Creation routine with Anisotropy
  /*!
   * \param fs 	    gauge state     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   * \param aniso   anisotropy struct   	       (Read)
   */
  void EvenOddPrecWilsonLinOp::create(Handle< FermState<T,P,Q> > fs,
				      const Real& Mass_,
				      const AnisoParam_t& anisoParam)
  {
    START_CODE();

    create(fs, Mass_, makeFermCoeffs(anisoParam));

    END_CODE();
  }


  //! Creation routine with general coefficients
  /*!
   * \param fs 	    gauge state     	       (Read)
   * \param Mass_   fermion kappa   	       (Read)
   * \param coeffs_ fermion coeffs 	       (Read)
   */
  void EvenOddPrecWilsonLinOp::create(Handle< FermState<T,P,Q> > fs,
				      const Real& Mass_,
				      const multi1d<Real>& coeffs_)
  {
    START_CODE();

    Mass   = Mass_;
    coeffs = coeffs_;

    if (coeffs.size() != Nd)
    {
      QDPIO::cerr << "EvenOddPrecWilsonLinOp::create : coeffs not size Nd" << endl;
      QDP_abort(1);
    }

    D.create(fs,coeffs);

    // Fold in the diagonal parts of the laplacian. Here, we insist that
    // the laplacian has the same scaling as the first deriv. term.
    fact = Mass;
    for(int mu=0; mu < Nd; ++mu)
      fact += coeffs[mu];

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

#if 0
 //! Derivative of even-odd linop component
  void 
  EvenOddPrecWilsonLinOp::derivEvenOddLinOpMP(multi1d<LatticeColorMatrix>& ds_u,
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
  EvenOddPrecWilsonLinOp::derivOddEvenLinOpMP(multi1d<LatticeColorMatrix>& ds_u,
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
#endif

  //! Return flops performed by the operator()
  unsigned long EvenOddPrecWilsonLinOp::nFlops() const
  { 
    unsigned long cbsite_flops = 2*D.nFlops()+6*Nc*Ns;
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

} // End Namespace Chroma
