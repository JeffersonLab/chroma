// -*- C++ -*-
// $Id: poly_cheb_fermact_w.h,v 3.3 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief Chebyshev polynomial fermion action
 */

#ifndef __poly_cheb_fermact_w_h__
#define __poly_cheb_fermact_w_h__

#include "wilstype_polyfermact_w.h"
#include "actions/ferm/linop/lgherm_w.h"

namespace Chroma
{
  //! Name and registration
  /*! \ingroup fermacts */
  namespace PolyChebFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Params for Chebyshev polynomial preconditioner
  /*! \ingroup fermacts */
  struct PolyChebFermActParams
  {
    PolyChebFermActParams() {}
    PolyChebFermActParams(XMLReader& in, const std::string& path);
    
    struct PolyParams
    {
      int degree;
      Real UpperBound ;
      Real LowerBound ;
      int order ;
    } polyParams;

    std::string  AuxFermAct;      /*!< string holding fermact xml */
  };

  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const string& path, PolyChebFermActParams& param);
  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const string& path, const PolyChebFermActParams& param);


  //! Chebyshev Polynomial fermion action
  /*! \ingroup fermacts
   *
   */
  class PolyChebFermAct : public PolyWilsonTypeFermAct<LatticeFermion, 
			  multi1d<LatticeColorMatrix> ,
			  multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    PolyChebFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
		    const PolyChebFermActParams& param_);

    //! Produce a linear operator for this action
    DiffLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const
      {
	return fermact->linOp(state);
      }

    //! Produce a linear operator M^dag.M for this action
    DiffLinearOperator<T,P,Q>* lMdagM(Handle< FermState<T,P,Q> > state) const
      {
	return fermact->lMdagM(state);
      }

    //! Produce the gamma_5 hermitian operator H_w
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const
      {
	return fermact->hermitianLinOp(state);
      }

    //! Produce a linear operator for this action
    DiffLinearOperator<T,P,Q>* polyPrecLinOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce a linear operator M^dag.M for this action
    PolyLinearOperator<T,P,Q>* polyLinOp(Handle< FermState<T,P,Q> > state) const;

    //! Return a linear operator solver for this action to solve M*psi=chi 
    /*! Default implementation provided */
    PolyPrecSystemSolver<T>* invPolyPrec(Handle< FermState<T,P,Q> > state,
					 const GroupXML_t& invParam) const;

    //! Destructor is automatic
    ~PolyChebFermAct() {}

  protected:
    //! Return the factory object that produces a state
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    /*! The user will supply the FermState in a derived class */
    PolyChebFermAct() {} //hide default constructor
    //! Assignment
    void operator=(const PolyChebFermAct& a) {}
   
  private:
    Handle< CreateFermState<T,P,Q> >  cfs;   /*!< fermion state creator */
    PolyChebFermActParams param;

    // A handle for the PrecWilsonFermAct
    Handle< WilsonTypeFermAct<T,P,Q> > fermact;
  };

}

#endif
