/*! \file
 *  \brief Even-odd preconditioned clover linear operator
 */

#include "actions/ferm/linop/eoprec_clover_linop_w.h"
#include "sse_dslash.h"

using namespace QDP::Hints;

namespace Chroma 
{ 

  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  fermion kappa   	       (Read)
   */
  void EvenOddPrecCloverLinOp::create(Handle< FermState<T,P,Q> > fs, 
				      const CloverFermActParams& param_)
  {
    START_CODE();

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << std::endl;

    param = param_;

    clov.create(fs, param);
 
    invclov.create(fs,param,clov);  // make a copy
    invclov.choles(0);  // invert the cb=0 part

    D.create(fs, param.anisoParam);

    clov_deriv_time = 0;
    clov_apply_time = 0;

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << std::endl;
    END_CODE();
  }

  //! Apply the the odd-odd block onto a source std::vector
  void 
  EvenOddPrecCloverLinOp::oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
				      enum PlusMinus isign) const
  {
    START_CODE();

    swatch.reset(); swatch.start();
    clov.apply(chi, psi, isign, 1);
    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();

    END_CODE();
  }


  //! Apply the the even-even block onto a source std::vector
  void 
  EvenOddPrecCloverLinOp::evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					enum PlusMinus isign) const
  {
    START_CODE();

    // Nuke for testing
    swatch.reset(); swatch.start();
    clov.apply(chi, psi, isign, 0);
    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();
    
    END_CODE();
  }

  //! Apply the inverse of the even-even block onto a source std::vector
  void 
  EvenOddPrecCloverLinOp::evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					   enum PlusMinus isign) const
  {
    START_CODE();

    swatch.reset(); swatch.start();
    invclov.apply(chi, psi, isign, 0);
    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();
    
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
  EvenOddPrecCloverLinOp::evenOddLinOp(LatticeFermion& chi, 
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
  EvenOddPrecCloverLinOp::oddEvenLinOp(LatticeFermion& chi, 
				       const LatticeFermion& psi, 
				       enum PlusMinus isign) const
  {
    START_CODE();

    Real mhalf = -0.5;

    D.apply(chi, psi, isign, 1);
    chi[rb[1]] *= mhalf;
  
    END_CODE();
  }


  //! Apply even-odd preconditioned Clover fermion linear operator
  /*!
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void EvenOddPrecCloverLinOp::operator()(LatticeFermion & chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
  {
    START_CODE();

    LatticeFermion tmp1; moveToFastMemoryHint(tmp1);
    LatticeFermion tmp2; moveToFastMemoryHint(tmp2);
    Real mquarter = -0.25;

    // Prepost receives for Dslash
    sse_su3dslash_prepost_receives();

    //  chi_o  =  A_oo  psi_o  -  tmp1_o
    clov.apply(chi, psi, isign, 1);
  
    //  tmp1_o  =  D_oe   A^(-1)_ee  D_eo  psi_o
    D.apply(tmp1, psi, isign, 0);

    // Prepost receives for Dslash
    sse_su3dslash_prepost_receives();
    invclov.apply(tmp2, tmp1, isign, 0);
    D.apply(tmp1, tmp2, isign, 1);

    //  chi_o  =  A_oo  psi_o  -  tmp1_o
    chi[rb[1]] += mquarter*tmp1;
    
    END_CODE();
  }


  //! Apply the even-even block onto a source std::vector
  void 
  EvenOddPrecCloverLinOp::derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					     const LatticeFermion& chi, const LatticeFermion& psi, 
					     enum PlusMinus isign) const
  {
    START_CODE();
    
    swatch.reset(); swatch.start();
    clov.deriv(ds_u, chi, psi, isign, 0);
    swatch.stop();
    clov_deriv_time  += swatch.getTimeInSeconds();

    END_CODE();
  }

  //! Apply the even-even block onto a source std::vector
  void 
  EvenOddPrecCloverLinOp::derivLogDetEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
						   enum PlusMinus isign) const
  {
    START_CODE();

    // Testing Odd Odd Term - get nothing from even even term
    invclov.derivTrLn(ds_u, isign, 0);
    
    END_CODE();
  }

  //! Apply the the even-odd block onto a source std::vector
  void 
  EvenOddPrecCloverLinOp::derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					    const LatticeFermion& chi, const LatticeFermion& psi, 
					    enum PlusMinus isign) const
  {
    START_CODE();
    ds_u.resize(Nd);
    D.deriv(ds_u, chi, psi, isign, 0);
    for(int mu=0; mu < Nd; mu++) { 
      ds_u[mu]  *= Real(-0.5);
    }
    END_CODE();
  }
 
  //! Apply the the odd-even block onto a source std::vector
  void 
  EvenOddPrecCloverLinOp::derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					    const LatticeFermion& chi, const LatticeFermion& psi, 
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

  // Inherit this
  //! Apply the the odd-odd block onto a source std::vector
  void 
  EvenOddPrecCloverLinOp::derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					   const LatticeFermion& chi, const LatticeFermion& psi, 
					   enum PlusMinus isign) const
  {   
    START_CODE();

    swatch.reset(); swatch.start();
    clov.deriv(ds_u, chi, psi, isign, 1);
    swatch.stop();
    clov_deriv_time += swatch.getTimeInSeconds();
    
    END_CODE();
  }


  //! Return flops performed by the operator()
  unsigned long EvenOddPrecCloverLinOp::nFlops() const
  {
    unsigned long cbsite_flops = 2*D.nFlops()+2*clov.nFlops()+4*Nc*Ns;
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

  //! Get the log det of the even even part
  // BUt for now, return zero for testing.
  Double EvenOddPrecCloverLinOp::logDetEvenEvenLinOp(void) const  {
    return invclov.cholesDet(0);
  }
} // End Namespace Chroma
