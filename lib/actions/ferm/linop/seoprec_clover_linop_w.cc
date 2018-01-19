/*! \file
 *  \brief Symmetric even-odd preconditioned clover linear operator
 */

#include "actions/ferm/linop/seoprec_clover_linop_w.h"

namespace Chroma 
{ 
  using namespace QDP::Hints;

  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  fermion kappa   	       (Read)
   */
  void SymEvenOddPrecCloverLinOp::create(Handle< FermState<T,P,Q> > fs, 
					 const CloverFermActParams& param_)
  {
    START_CODE();
    // QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << std::endl;

    param = param_;

    clov.create(fs, param);
 
    invclov.create(fs,param,clov);  // make a copy
    invclov.choles(0);  // invert the cb=0 part
    invclov.choles(1);  // invert the cb=1 part

    D.create(fs, param.anisoParam);

    clov_deriv_time = 0;
    clov_apply_time = 0;

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << std::endl;
    END_CODE();
  }

  //! Apply the the odd-odd block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					 enum PlusMinus isign) const
  {
    START_CODE();

    swatch.reset(); swatch.start();
    clov.apply(chi, psi, isign, 1);
    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();

    END_CODE();
  }


  //! Apply the inverse of the odd-odd block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::oddOddInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					    enum PlusMinus isign) const
  {
    START_CODE();

    swatch.reset(); swatch.start();
    invclov.apply(chi, psi, isign, 1);
    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();
    
    END_CODE();
  }
  

  //! Apply the the even-even block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
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
  SymEvenOddPrecCloverLinOp::evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
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
  SymEvenOddPrecCloverLinOp::evenOddLinOp(LatticeFermion& chi, 
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
  SymEvenOddPrecCloverLinOp::oddEvenLinOp(LatticeFermion& chi, 
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
  void SymEvenOddPrecCloverLinOp::operator()(LatticeFermion & chi, 
					     const LatticeFermion& psi, 
					     enum PlusMinus isign) const
  {
    START_CODE();

    LatticeFermion tmp1; moveToFastMemoryHint(tmp1);
    LatticeFermion tmp2; moveToFastMemoryHint(tmp2);
    Real mquarter = -0.25;
  
    //  tmp2_o  =  A^(-1)_oo  D_oe  A^(-1)_ee  D_eo  psi_o
    D.apply(tmp1, psi, isign, 0);

    swatch.reset(); swatch.start();
    invclov.apply(tmp2, tmp1, isign, 0);
    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();

    D.apply(tmp1, tmp2, isign, 1);

    swatch.reset(); swatch.start();
    invclov.apply(tmp2, tmp1, isign, 1);
    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();

    //  chi_o  =  psi_o  -  tmp2_o
    chi[rb[1]] = psi  +  mquarter*tmp2;

    // Twisted Term?
    if( param.twisted_m_usedP ){ 
      QDPIO::cerr << __func__ << ": do not support twisted mass terms" << std::endl;
      QDP_abort(1);
    }

    END_CODE();
  }


  //! Apply the even-even block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
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
  SymEvenOddPrecCloverLinOp::derivLogDetEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
						      enum PlusMinus isign) const
  {
    START_CODE();

    invclov.derivTrLn(ds_u, isign, 0);
    
    END_CODE();
  }

  //! Apply the even-even block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::derivLogDetOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u,
						    enum PlusMinus isign) const
  {
    START_CODE();

    invclov.derivTrLn(ds_u, isign, 1);
    
    END_CODE();
  }

  //! Apply the the even-odd block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
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
  SymEvenOddPrecCloverLinOp::derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
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
  SymEvenOddPrecCloverLinOp::derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
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
  unsigned long SymEvenOddPrecCloverLinOp::nFlops() const
  {
    unsigned long cbsite_flops = 2*D.nFlops()+2*clov.nFlops()+4*Nc*Ns;
    if(  param.twisted_m_usedP ) { 
      cbsite_flops += 4*Nc*Ns; // a + mu*b : a = chi, b = g_5 I psi
    }
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

  //! Get the log det of the even even part
  // BUt for now, return zero for testing.
  Double SymEvenOddPrecCloverLinOp::logDetEvenEvenLinOp(void) const  {
    return invclov.cholesDet(0);
  }

  //! Get the log det of the odd odd part
  // BUt for now, return zero for testing.
  Double SymEvenOddPrecCloverLinOp::logDetOddOddLinOp(void) const  {
    return invclov.cholesDet(1);
  }

} // End Namespace Chroma
