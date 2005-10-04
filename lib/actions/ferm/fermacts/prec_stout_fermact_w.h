// -*- C++ -*-
// $Id: prec_stout_fermact_w.h,v 2.1 2005-10-04 19:23:19 bjoo Exp $

/*! @file
 *  @brief Proxy fermion action class instance for unpreconditioned stout fermacts 
 */

#ifndef _prec_stout_fermact_w_h_
#define _prec_stout_fermact_w_h_

#include "chromabase.h"
#include "fermact.h"
#include "actions/ferm/fermacts/stout_state.h"
#include "actions/ferm/fermacts/stout_fermact_params.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"

namespace Chroma 
{


  //! Name and registration
  namespace EvenOddPrecStoutWilsonTypeFermActEnv {
    extern const std::string name;
    extern const bool registered;
  }


  //! A proxy class for unpreconditioned stout actions
  /*! @ingroup actions 
   *
   * Basically this shell of an action substitutes a call
   * to create state to create a stout state
   *
   * apart from this it does only gymnastics to forward 
   * to an embedded fermion action 
   */
  
  class EvenOddPrecStoutWilsonTypeFermAct 
    : public EvenOddPrecWilsonTypeFermAct<LatticeFermion , 
                                     multi1d<LatticeColorMatrix> > 
    {

    public:
      // Destructor is automagic
      ~EvenOddPrecStoutWilsonTypeFermAct() {}

      //! General FermBC
      EvenOddPrecStoutWilsonTypeFermAct(Handle< FermBC<LatticeFermion> > fbc_,
			 const StoutFermActParams& p_):  p(p_) { 

	struct EvenOddPrecCastFailure {
	  EvenOddPrecCastFailure(std::string e) : auxfermact(e) {};
	  const string auxfermact;
	};

	try {
	  // Make an XML Reader from the internal fermact
	  std::istringstream is(p_.internal_fermact);
	  XMLReader fermact_xml(is);
	  
	  // Get to the top of it
	  XMLReader internal_reader(fermact_xml, "/InternalFermionAction");
	  
	  // Get the name of the internal fermion action
	  std::string if_name;
	  read(internal_reader, "FermAct", if_name);
	  
	  // Creaate it using the FermionActionFactory
	  FermionAction<LatticeFermion>* FA_ptr= TheFermionActionFactory::Instance().createObject(if_name, fermact_xml, "/InternalFermionAction"); 

	  // Upcast to EvenOddPrecWilsonTypeFermact
	  EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >* S_internal = dynamic_cast< EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >* >(FA_ptr);

	  
	  if( S_internal == 0x0 ) { 
	    throw EvenOddPrecCastFailure(p_.internal_fermact);
	  }
	  else {
	    S_w = S_internal;
	  }
	}
	catch(const std::string& e) {
	  QDPIO::cout << "Error creating Internal Fermact: " << e << endl;
	  QDP_abort(1);
	}
	catch(const EvenOddPrecCastFailure& e) {
	  QDPIO::cout << "Internal Fermact cannot be cast to EvenOddPrecWilsonTypeFermAct: " << e.auxfermact << endl;
	  QDP_abort(1);
	}

      }

      //! Copy Constructor 
      EvenOddPrecStoutWilsonTypeFermAct(const EvenOddPrecStoutWilsonTypeFermAct& a) :
							p(a.p),
							S_w(a.S_w) { }

      //! Assignment
      EvenOddPrecStoutWilsonTypeFermAct& operator=(const EvenOddPrecStoutWilsonTypeFermAct& a) 
      {
	p = a.p;
	S_w = a.S_w;
	return *this;
      }

      const FermBC<LatticeFermion>& getFermBC() const { 
	return S_w->getFermBC();
      }

      //! Produce a linear operator for this action
      const EvenOddPrecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const 
      {	
	return S_w->linOp(state);
      }


      //! Produce a linear operator M^dag.M for this action
      const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const {
	return S_w->lMdagM(state);
      }

      //! Produce the gamma_5 hermitian operator H_w
      const LinearOperator<LatticeFermion>* hermitianLinOp(Handle< const ConnectState> state) const { 
	return S_w->hermitianLinOp(state);
      }


      const ConnectState* createState(const multi1d<LatticeColorMatrix>& u_thin) const
      {
	multi1d<LatticeColorMatrix> u_tmp(Nd);
	u_tmp = u_thin;

	// Apply BC-s to thin links
	S_w->getFermBC().modifyU(u_tmp);

	// Make a stout state
	return new StoutConnectState(u_tmp, p.rho, p.n_smear);
      }

      const SystemSolver<LatticeFermion>* qprop(Handle<const ConnectState> state, 
						const InvertParam_t& invParam) const {
	return S_w->qprop(state, invParam);
      }

    private:
      StoutFermActParams p;

      Handle< EvenOddPrecWilsonTypeFermAct<LatticeFermion, 
	                              multi1d<LatticeColorMatrix> > > S_w;

  }; 
};


#endif
