// -*- C++ -*-
// $Id: eoprec_twm_fermact_array_w.h,v 1.1 2008-11-04 18:42:58 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Twisted-mass where each flavor is one of two array elements
 */

#ifndef __prec_twm_fermact_array_w_h__
#define __prec_twm_fermact_array_w_h__

#include "eoprec_constdet_wilstype_fermact_w.h"

namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecTwmFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Params for wilson ferm acts
  struct EvenOddPrecTwmFermActArrayParams
  {
    EvenOddPrecTwmFermActArrayParams() {}
    EvenOddPrecTwmFermActArrayParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    Real mu_sigma;
    Real mu_delta;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, EvenOddPrecTwmFermActArrayParams& param);
  void write(XMLWriter& xml, const string& path, const EvenOddPrecTwmFermActArrayParams& param);


  //! Even-odd preconditioned Wilson fermion action with parity breaking term
  /*! \ingroup fermacts
   *
   * Even-odd preconditioned wilson fermion action with parity breaking term
   * Only defined on odd subset.
   *
   * The kernel for Wilson fermions with a parity breaking term is
   *
   *      M  =  (d+M) + i*mu_sigma*gamma_5*tau_1 + mu_delta*tau_3  - (1/2) D'
   */
  class EvenOddPrecTwmFermActArray : public EvenOddPrecConstDetWilsonTypeFermAct5D<LatticeFermion, 
				     multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    EvenOddPrecTwmFermActArray(Handle< CreateFermState<T,P,Q> > fs_, 
			       const EvenOddPrecTwmFermActArrayParams& param_) : 
      fs(fs_), param(param_) {}

    //! Copy constructor
    EvenOddPrecTwmFermActArray(const EvenOddPrecTwmFermActArray& a) : 
      fs(a.fs), param(a.param) {}

    //! Length of DW flavor index/space
    int size() const {return 2;}

    //! Produce a linear operator for this action
    EvenOddPrecConstDetLinearOperatorArray<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Override to produce an even-odd prec. Pauli-Villars linear operator for this action
    EvenOddPrecConstDetLinearOperatorArray<T,P,Q>* linOpPV(Handle< FermState<T,P,Q> > state) const
      {
	QDP_error_exit("linOpPV not implemented yet for this action\n");
	return 0;
      }

    //! Produce the gamma_5 hermitin op gamma_5 M
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      {
	QDP_error_exit("gamma5HermLinOp not implemented yet for this action\n");
	return 0;
      }

    //! Return quark prop solver, solution of unpreconditioned system
    SystemSolverArray<T>* qpropT(Handle< FermState<T,P,Q> > state,
				 const GroupXML_t& invParam) const;

    //! Destructor is automatic
    ~EvenOddPrecTwmFermActArray() {}

  protected:
    //! Return the fermion create state for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *fs;}

    //! Partial constructor
    EvenOddPrecTwmFermActArray() {}
    //! Hide =
    void operator=(const EvenOddPrecTwmFermActArray& a) {}

  private:
    Handle< CreateFermState<T,P,Q> >  fs;
    EvenOddPrecTwmFermActArrayParams param;
  };

}


#endif
