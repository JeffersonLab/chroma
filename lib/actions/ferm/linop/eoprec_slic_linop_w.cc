// $Id: eoprec_slic_linop_w.cc,v 3.4 2009-05-19 20:07:33 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned clover linear operator
 */

#include "actions/ferm/linop/eoprec_slic_linop_w.h"

using namespace QDP::Hints;

/*
           Notes on SLIC

The idea is you want *two* optimized dslash calls of the form

dslash(0.5*U + 0.5*U',PLUS)*psi
dslash(0.5*U - 0.5*U',MINUS)*psi

If you work this out, you'll find the cross terms cancel, and
you are left with

Dslash' psi(x) = sum_mu  [ gamma_mu U_mu(x) - U'_mu(x) ] psi(x + mu)
                      - [ gamma_mu U^\dag_mu(x - mu) + U'^\dag_mu(x - mu) ] psi(x-mu)

The point is the "MINUS" in the dslash call (the dagger) flips the
sign and changes  gamma_mu*U - U'   ->  gamma_mu*U + U'

So, two copies of the gauge fields are needed beyond
holding the fat and the thin links, namely

U_1 = 0.5*U_thin + 0.5*U_fat
U_2 = 0.5*U_thin - 0.5*U_fat

and then one does (in pseudo-code)

chi  = dslash(U_1,PLUS)*psi
chi += dslash(U_2,MINUS)*psi

So, to reiterate. The existing optimized codes can be used. You
get the projector trick. You need 2 dslash calls, so SLIC (or SLRC)
is twice the cost of normal Wilson ferms that are fully fattened
or not fattened at all. However, this factor of two is better than
the naive cost of *four* which would come from doing a full
3x3(gauge) * 3x4(full spinor) . The spin projector trick allows us
to reduce the actual work to
3x3(gauge) * 3x2(half spinor)

This factor of 2 overhead is what Waseem advertises as the normal
cost of SLIC.
*/


namespace Chroma 
{ 

  //! Full constructor
  EvenOddPrecSLICLinOp::EvenOddPrecSLICLinOp(Handle< FermState<T,P,Q> > fs,
					     const CloverFermActParams& param_) : slic_fs(fs)
  {
    create(fs, param_);
  }



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

    multi1d<LatticeColorMatrix> u_pl(Nd);
    multi1d<LatticeColorMatrix> u_mn(Nd);

    for(int mu=0; mu < Nd; ++mu)
    {
      u_pl[mu] = (fs.cast<SLICFermState<T,P,Q> >())->getThinLinks()[mu]
	+ (fs.cast<SLICFermState<T,P,Q> >())->getLinks()[mu];

      u_mn[mu] = (fs.cast<SLICFermState<T,P,Q> >())->getThinLinks()[mu]
	- (fs.cast<SLICFermState<T,P,Q> >())->getLinks()[mu];

      u_pl[mu] *= Real(0.5);
      u_mn[mu] *= Real(0.5);
    }

    // Need to make sure that fs is a stout ferm state
    clov.create(fs, param);
 
    invclov.create(fs,param,clov);  // make a copy
    invclov.choles(0);  // invert the cb=0 part

    // Create the two dslashs
    //
    // WARNING: there are factor of 1/2 floating around. Need to fix this in the
    //          derivatives below.
    // 
    fs_pl  = new SimpleFermState<T,P,Q>(fs->getFermBC(), u_pl);
    fs_mn  = new SimpleFermState<T,P,Q>(fs->getFermBC(), u_mn);

    D_pl.create(fs_pl, param.anisoParam);
    D_mn.create(fs_mn, param.anisoParam);
    
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

    enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;
    Real mhalf = -0.5;

    LatticeFermion tmp;
    D_pl.apply(chi, psi, isign, 0);
    D_mn.apply(tmp, psi, msign, 0);
    chi[rb[0]] += tmp;
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

    enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;
    Real mhalf = -0.5;

    LatticeFermion tmp;
    D_pl.apply(chi, psi, isign, 1);
    D_mn.apply(tmp, psi, msign, 1);
    chi[rb[1]] += tmp;
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

    enum PlusMinus msign = (isign == PLUS) ? MINUS : PLUS;

    LatticeFermion tmp1; moveToFastMemoryHint(tmp1);
    LatticeFermion tmp2; moveToFastMemoryHint(tmp2);
    LatticeFermion tmp3; moveToFastMemoryHint(tmp3);
    Real mquarter = -0.25;
  
    //  tmp1_o  =  D_oe   A^(-1)_ee  D_eo  psi_o
    D_pl.apply(tmp1, psi, isign, 0);
    D_mn.apply(tmp3, psi, msign, 0);
    tmp1[rb[0]] += tmp3;
    invclov.apply(tmp2, tmp1, isign, 0);
    D_pl.apply(tmp1, tmp2, isign, 1);
    D_mn.apply(tmp3, tmp2, msign, 1);
    tmp1[rb[1]] += tmp3;

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

    //
    //
    // WARNING: these derivatives are not correct
    //
    //
    std::cerr << __func__ << ": need to fix up SLIC derivatives\n";
    QDP_abort(1);


    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    ds_u.resize(Nd);

    clov.deriv(ds_tmp, chi, psi, isign, 0);
    SLICFermState<T,P,Q>& sfs = dynamic_cast<SLICFermState<T,P,Q>& >(*slic_fs);

    sfs.fatForceToThin(ds_tmp,ds_u);
    
    END_CODE();
  }

  //! Apply the even-even block onto a source vector
  void 
  EvenOddPrecSLICLinOp::derivLogDetEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
						 enum PlusMinus isign) const
  {
    START_CODE();
    
    // Testing Odd Odd Term - get nothing from even even term
    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    ds_u.resize(Nd);

    invclov.derivTrLn(ds_tmp, isign, 0);
    SLICFermState<T,P,Q>& sfs = dynamic_cast<SLICFermState<T,P,Q>& >(*slic_fs);

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
    D_pl.deriv(ds_u, chi, psi, isign, 0);
//    D_mn.deriv(ds_u, chi, psi, isign, 0);

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

    D_pl.deriv(ds_u, chi, psi, isign, 1);
//    D_mn.deriv(ds_u, chi, psi, isign, 1);

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
    unsigned long cbsite_flops = 2*2*D_pl.nFlops()+2*clov.nFlops()+4*Nc*Ns;
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

  //! Get the log det of the even even part
  // BUt for now, return zero for testing.
  Double EvenOddPrecSLICLinOp::logDetEvenEvenLinOp(void) const  {
    return invclov.cholesDet(0);
  }
} // End Namespace Chroma
