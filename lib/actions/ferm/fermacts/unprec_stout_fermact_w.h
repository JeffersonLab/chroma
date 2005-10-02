// -*- C++ -*-
// $Id: unprec_stout_fermact_w.h,v 2.2 2005-10-02 19:33:12 bjoo Exp $

/*! @file
 *  @brief Proxy fermion action class instance for unpreconditioned stout fermacts 
 */

#ifndef _unprec_stout_fermact_w_h_
#define _unprec_stout_fermact_w_h_

#include "chromabase.h"
#include "fermact.h"
#include "actions/ferm/fermacts/stout_state.h"
#include "actions/ferm/fermacts/stout_fermact_params.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"

namespace Chroma 
{


  //! Name and registration
  namespace UnprecStoutWilsonTypeFermActEnv {
    extern const std::string name;
    extern const bool registered;
    WilsonTypeFermAct<LatticeFermion,multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,const std::string& path);
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
  
  class UnprecStoutWilsonTypeFermAct 
    : public UnprecWilsonTypeFermAct<LatticeFermion , 
                                     multi1d<LatticeColorMatrix> > 
    {

    public:
      // Destructor is automagic
      ~UnprecStoutWilsonTypeFermAct() {}

      //! General FermBC
      UnprecStoutWilsonTypeFermAct(Handle< FermBC<LatticeFermion> > fbc_,
			 const StoutFermActParams& p_):  p(p_) { 

	struct UnprecCastFailure {
	  UnprecCastFailure(std::string e) : auxfermact(e) {};
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

	  // Upcast to UnprecWilsonTypeFermact
	  UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >* S_internal = dynamic_cast< UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >* >(FA_ptr);

	  
	  if( S_internal == 0x0 ) { 
	    throw UnprecCastFailure(p_.internal_fermact);
	  }
	  else {
	    S_w = S_internal;
	  }
	}
	catch(const std::string& e) {
	  QDPIO::cout << "Error creating Internal Fermact: " << e << endl;
	  QDP_abort(1);
	}
	catch(const UnprecCastFailure& e) {
	  QDPIO::cout << "Internal Fermact cannot be cast to UnprecWilsonTypeFermAct: " << e.auxfermact << endl;
	  QDP_abort(1);
	}

      }

      //! Copy Constructor 
      UnprecStoutWilsonTypeFermAct(const UnprecStoutWilsonTypeFermAct& a) :
							p(a.p),
							S_w(a.S_w) { }

      //! Assignment
      UnprecStoutWilsonTypeFermAct& operator=(const UnprecStoutWilsonTypeFermAct& a) 
      {
	p = a.p;
	S_w = a.S_w;
	return *this;
      }

      const FermBC<LatticeFermion>& getFermBC() const { 
	return S_w->getFermBC();
      }

      //! Produce a linear operator for this action
      const UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const 
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

    private:
      StoutFermActParams p;

      Handle< UnprecWilsonTypeFermAct<LatticeFermion, 
	                              multi1d<LatticeColorMatrix> > > S_w;

  };
};


#endif
