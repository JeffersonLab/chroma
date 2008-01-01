// -*- C++ -*-
// $Id: unprec_graphene_linop_w.h,v 1.2 2008-01-01 22:13:06 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Graphene fermion linear operator.
 *
 * This formulation follows Borici's variant of Creutz's graphene
 * fermion construction. Borici's variant is described in
 * arXiv:0712.4401 and Cruetz's original construction is described
 * in arXiv:0712.1201
 */

#ifndef __unprec_graphene_linop_w_h__
#define __unprec_graphene_linop_w_h__

#include "linearop.h"
#include "state.h"
#include "io/aniso_io.h"

namespace Chroma 
{ 
  //! Unpreconditioned Graphene operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   *                                                      ~      ~+
   * This subroutine applies the unpreconditioned matrix  M  or  M   the vector
   * Psi,
   *
   *      	       	   {   ~
   *      	       	   {   M(U) . Psi      	       if  ISign = PLUS
   *      	   Chi  =  {
   *      	       	   {   ~   +
   *      	       	   {   M(U)  . Psi     	       if  ISign = MINUS
   *
   * Algorithm:
   *
   * The kernel for these Graphene fermions is
   *
   *      M  =  m + i\sum_\mu\gamma_mu + 
   *          (1/2)*[(i*\Gamma_mu+\gamma_mu)*H_\mu + (i*\Gamma_mu-\gamma_mu)*H_-\mu]
   *
   */
  
  class UnprecGrapheneLinOp : public UnprecLinearOperator<LatticeFermion, 
                     multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    UnprecGrapheneLinOp() {}

    //! Full constructor
    UnprecGrapheneLinOp(Handle< FermState<T,P,Q> > fs, const Real& Mass_)
      {create(fs,Mass_);}

    //! Full constructor with Anisotropy
    UnprecGrapheneLinOp(Handle< FermState<T,P,Q> > fs,
			const Real& Mass_,
			const AnisoParam_t& aniso)
      {create(fs,Mass_,aniso);}

    //! Destructor is automatic
    ~UnprecGrapheneLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return *fbc;}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& Mass_);

    //! Creation routine with Anisotropy
    void create(Handle< FermState<T,P,Q> > fs,
		const Real& Mass_,
		const AnisoParam_t& aniso);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

    //! Derivative of unpreconditioned Graphene dM/dU
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const LatticeFermion& chi, const LatticeFermion& psi, 
	       enum PlusMinus isign) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  private:
    //! Form  gamma_mu * psi
    void gammaMults(multi1d<LatticeFermion>& tmp1, const LatticeFermion& psi) const;

    //! Form  i*Gamma_mu * psi
    void iGamMu(LatticeFermion& iGam, const multi1d<LatticeFermion>& gams, int mu) const;

  private:
    Real             Mass;
    AnisoParam_t     anisoParam;

    Handle< FermBC<T,P,Q> >       fbc;
    multi1d<LatticeColorMatrix>   u;
    multi2d<int>                  alpha;
  };

} // End Namespace Chroma


#endif
