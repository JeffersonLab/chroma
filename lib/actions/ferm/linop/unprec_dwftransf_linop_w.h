// -*- C++ -*-
// $Id: unprec_dwftransf_linop_w.h,v 1.1 2004-11-02 10:33:50 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion linear operator
 */

#ifndef __unprec_dwftransf_linop_w_h__
#define __unprec_dwftransf_linop_w_h__

#include "linearop.h"
#include "handle.h"
#include "invtype.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"

using namespace QDP;
using namespace Chroma;

class UnprecDWFTransfLinOp : public LinearOperator<LatticeFermion>
{
public:
  //! Partial constructor
  UnprecDWFTransfLinOp() {}

  //! Full constructor
  UnprecDWFTransfLinOp(const multi1d<LatticeColorMatrix>& u_, 
		      const Real& Mass_,
		      const Real& b5_,
		      const Real& c5_,
		      const InvertParam_t& invParam_)
    {create(u_,Mass_, b5_, c5_, invParam_);}

  //! Destructor is automatic
  ~UnprecDWFTransfLinOp() {}

  //! Only defined on the odd subset
  const OrderedSubset& subset() const {return all;}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u_, 
	      const Real& Mass_,
	      const Real& b5_,
	      const Real& c5_,
	      const InvertParam_t& invParam_);

  //! Apply the operator onto a source vector
  void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;


private:
  Real Mass;
  Real b5;
  Real c5;
  InvertParam_t inv_param;
  multi1d<LatticeColorMatrix> u;
  Handle< LinearOperator<LatticeFermion> > D_w;
  Handle< LinearOperator<LatticeFermion> > D_denum;  
  Handle< LinearOperator<LatticeFermion> > D_sq_denum;
};


class UnprecDWFTransfMdagMLinOp : public LinearOperator<LatticeFermion>
{
public:
  //! Partial constructor
  UnprecDWFTransfMdagMLinOp() {}

  //! Full constructor
  UnprecDWFTransfMdagMLinOp(const multi1d<LatticeColorMatrix>& u_, 
			    const Real& Mass_,
			    const Real& b5_,
			    const Real& c5_,
			    const InvertParam_t& invParam_)
  {create(u_,Mass_, b5_, c5_, invParam_);}

  //! Destructor is automatic
  ~UnprecDWFTransfMdagMLinOp() {}

  //! Only defined on the odd subset
  const OrderedSubset& subset() const {return all;}

  //! Creation routine
  void create(const multi1d<LatticeColorMatrix>& u_, 
	      const Real& Mass_,
	      const Real& b5_,
	      const Real& c5_,
	      const InvertParam_t& invParam_);

  //! Apply the operator onto a source vector
  void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

private:
  Real Mass;
  Real b5;
  Real c5;
  InvertParam_t inv_param;
  multi1d<LatticeColorMatrix> u;
  Handle< LinearOperator<LatticeFermion> > D_w;
  Handle< LinearOperator<LatticeFermion> > D_sq_denum;

};


// Operators to apply the denominator
class UnprecDWFTransfDenLinOp : public LinearOperator<LatticeFermion>
{
public:
  //! Partial constructor
  UnprecDWFTransfDenLinOp() {}

  //! Full constructor
  UnprecDWFTransfDenLinOp(const multi1d<LatticeColorMatrix>& u_, 
			 const Real& b5_minus_c5_,
			 const Handle<LinearOperator<LatticeFermion> > D_w_ ) 
    : u(u_), b5_minus_c5(b5_minus_c5_), D_w(D_w_) {}


  //! Destructor is automatic
  ~UnprecDWFTransfDenLinOp() {}

  //! Only defined on the odd subset
  const OrderedSubset& subset() const {return all;}

  //! Apply the operator onto a source vector
  void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const 
  {

    // Apply chi ( 2 + (b5-c5) D_w ) psi
    (*D_w)(chi, psi, isign);
    chi *= b5_minus_c5;
    chi += Real(2)*psi;
  }


private:
  multi1d<LatticeColorMatrix> u;
  Real b5_minus_c5;
  Handle< LinearOperator<LatticeFermion> > D_w;
};

class UnprecDWFTransfMdagMDenLinOp : public LinearOperator<LatticeFermion>
{
public:
  //! Partial constructor
  UnprecDWFTransfMdagMDenLinOp() {}

  //! Full constructor
  UnprecDWFTransfMdagMDenLinOp(const multi1d<LatticeColorMatrix>& u_, 
			      const Real& b5_minus_c5_,
			      const Handle<LinearOperator<LatticeFermion> > D_w_ ) 
    : u(u_), b5_minus_c5(b5_minus_c5_), D_w(D_w_) {}
  

  //! Destructor is automatic
  ~UnprecDWFTransfMdagMDenLinOp() {}

  //! Only defined on the odd subset
  const OrderedSubset& subset() const {return all;}

  //! Apply the operator onto a source vector
  void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const 
  {

    // Apply chi (2 + (b5-c5)D^{dag})( 2 + (b5-c5) D_w ) psi

    // First apply tmp = (2 + (b5-c5)D_w) psi
    LatticeFermion tmp;
    (*D_w)(tmp, psi, isign);
    tmp *= b5_minus_c5;
    tmp += Real(2)*psi;

    // Now apply chi = ( 2 + (b5-c5)D_w^{dag}) tmp
    enum PlusMinus isign_dag = (isign == PLUS) ? MINUS : PLUS ;
    (*D_w)(chi, tmp, isign_dag);
    chi *= b5_minus_c5;
    chi +=Real(2)*tmp;
  }


private:
  multi1d<LatticeColorMatrix> u;
  Real b5_minus_c5;
  Handle< LinearOperator<LatticeFermion> > D_w;
};

#endif
