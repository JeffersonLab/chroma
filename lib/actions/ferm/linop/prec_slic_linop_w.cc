// $Id: prec_slic_linop_w.cc,v 1.3 2006-09-06 19:05:10 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned clover linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_slic_linop_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 

  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  fermion kappa   	       (Read)
   */
  void EvenOddPrecSLICLinOp::create(Handle< FermState<T,P,Q> > fs, 
				    const CloverFermActParams& param_)
  {
    START_CODE();

    param = param_;

     // Need to make sure that fs is a stout ferm state
    thin_fs  = new SimpleFermState<T,P,Q>( fs->getFermBC(), (fs.cast<SLICFermState<T, P, Q> >())->getThinLinks());
    clov.create(fs, param_);
 

    invclov = clov;  // make a copy
    invclov.choles(0);  // invert the cb=0 part

    D.create(thin_fs.cast< FermState<T,P,Q> >(), param_.anisoParam);
    
    END_CODE();
  }

  //! Apply the the odd-odd block onto a source vector
  void 
  EvenOddPrecSLICLinOp::oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
				      enum PlusMinus isign) const
  {
    clov.apply(chi, psi, isign, 1);
  }


  //! Apply the the even-even block onto a source vector
  void 
  EvenOddPrecSLICLinOp::evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					enum PlusMinus isign) const
  {
    // Nuke for testing
    clov.apply(chi, psi, isign, 0);
  }

  //! Apply the inverse of the even-even block onto a source vector
  void 
  EvenOddPrecSLICLinOp::evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					   enum PlusMinus isign) const
  {
    invclov.apply(chi, psi, isign, 0);
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
  EvenOddPrecSLICLinOp::evenOddLinOp(LatticeFermion& chi, 
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
  EvenOddPrecSLICLinOp::oddEvenLinOp(LatticeFermion& chi, 
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
  void EvenOddPrecSLICLinOp::operator()(LatticeFermion & chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
  {
    START_CODE();

    LatticeFermion tmp1; moveToFastMemoryHint(tmp1);
    LatticeFermion tmp2; moveToFastMemoryHint(tmp2);
    Real mquarter = -0.25;
  
    //  tmp1_o  =  D_oe   A^(-1)_ee  D_eo  psi_o
    D.apply(tmp1, psi, isign, 0);
    invclov.apply(tmp2, tmp1, isign, 0);
    D.apply(tmp1, tmp2, isign, 1);

    //  chi_o  =  A_oo  psi_o  -  tmp1_o
    clov.apply(chi, psi, isign, 1);

    chi[rb[1]] += mquarter*tmp1;
    
    END_CODE();
  }


  //! Apply the even-even block onto a source vector
  void 
  EvenOddPrecSLICLinOp::derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					     const LatticeFermion& chi, const LatticeFermion& psi, 
					     enum PlusMinus isign) const
  {
    START_CODE();

    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    ds_u.resize(Nd);

    clov.deriv(ds_tmp, chi, psi, isign, 0);
    SLICFermState<T, P, Q>& sfs = dynamic_cast<SLICFermState<T,P,Q>& >(*slic_fs);

    sfs.fatForceToThin(ds_tmp,ds_u);
    
    END_CODE();
  }

  //! Apply the even-even block onto a source vector
  void 
  EvenOddPrecSLICLinOp::derivEvenEvenLogDet(multi1d<LatticeColorMatrix>& ds_u,
					      enum PlusMinus isign) const
  {
    START_CODE();
    
    // Testing Odd Odd Term - get nothing from even even term
    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    ds_u.resize(Nd);

    invclov.derivTrLn(ds_tmp, isign, 0);
    SLICFermState<T, P, Q>& sfs = dynamic_cast<SLICFermState<T,P,Q>& >(*slic_fs);

    sfs.fatForceToThin(ds_tmp,ds_u);
    
    END_CODE();
  }

  //! Apply the the even-odd block onto a source vector
  void 
  EvenOddPrecSLICLinOp::derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					    const LatticeFermion& chi, const LatticeFermion& psi, 
					    enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);
    D.deriv(ds_u, chi, psi, isign, 0);

    // Undo ferm boundaries on ds_U
    slic_fs->getFermBC()->modify(ds_u);

    // But reinforce gauge boundaries
    slic_fs->getFermBC()->zero(ds_u);

    for(int mu=0; mu < Nd; mu++) { 
      ds_u[mu]  *= Real(-0.5);
    }
    END_CODE();
  }
 
  //! Apply the the odd-even block onto a source vector
  void 
  EvenOddPrecSLICLinOp::derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					    const LatticeFermion& chi, const LatticeFermion& psi, 
					    enum PlusMinus isign) const
  {
    START_CODE();
    ds_u.resize(Nd);

    D.deriv(ds_u, chi, psi, isign, 1);

    // Undo ferm boundaries on ds_U
    slic_fs->getFermBC()->modify(ds_u);

    // But reinforce gauge boundaries
    slic_fs->getFermBC()->zero(ds_u);

    for(int mu=0; mu < Nd; mu++) { 
     ds_u[mu]  *= Real(-0.5);
    }
    END_CODE();
  }

  // Inherit this
  //! Apply the the odd-odd block onto a source vector
  void 
  EvenOddPrecSLICLinOp::derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					   const LatticeFermion& chi, const LatticeFermion& psi, 
					   enum PlusMinus isign) const
  {   
    START_CODE();

    ds_u.resize(Nd);
    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    clov.deriv(ds_tmp, chi, psi, isign, 1);
    SLICFermState<T,P,Q>& sfs = dynamic_cast< SLICFermState<T,P,Q>& >(*slic_fs);
    sfs.fatForceToThin(ds_tmp,ds_u);
    
    END_CODE();
  }


  //! Return flops performed by the operator()
  unsigned long EvenOddPrecSLICLinOp::nFlops() const
  {
    unsigned long cbsite_flops = 2*D.nFlops()+2*clov.nFlops()+4*Nc*Ns;
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

  //! Get the log det of the even even part
  // BUt for now, return zero for testing.
  Double EvenOddPrecSLICLinOp::LogDetEvenEven(void) const  {
    return invclov.cholesDet(0);
  }
} // End Namespace Chroma
