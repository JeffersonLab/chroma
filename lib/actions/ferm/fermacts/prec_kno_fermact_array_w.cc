// $Id: prec_kno_fermact_array_w.cc,v 1.3 2004-11-23 20:02:50 kostas Exp $
/*! \file
 *  \brief preconditioned KNO fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_kno_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"
#include "actions/ferm/linop/prec_nef_general_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/fermfactory_w.h"
#include "actions/ferm/fermacts/zolotarev.h"

using namespace Chroma;

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecKNOFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct< multi1d<LatticeFermion> >* createFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
								XMLReader& xml_in,
								const std::string& path)
    {
      return new EvenOddPrecKNOFermActArray(fbc, EvenOddPrecKNOFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    EvenOddPrecDWFermActBaseArray<LatticeFermion>* createDWFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
							      XMLReader& xml_in,
							      const std::string& path)
    {
      return new EvenOddPrecKNOFermActArray(fbc, EvenOddPrecKNOFermActArrayParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "KNO";

    //! Register the Wilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActArrayFactory::Instance().registerObject(name, createFermAct)
                          & Chroma::TheEvenOddPrecDWFermActBaseArrayFactory::Instance().registerObject(name, createDWFermAct); 
  }


  //! Read parameters
  EvenOddPrecKNOFermActArrayParams::EvenOddPrecKNOFermActArrayParams(XMLReader& xml, 
							   const std::string& path)
  {
    XMLReader paramtop(xml, path);
    try {
      // Read the stuff for the action
      read(paramtop, "OverMass", OverMass);
      read(paramtop, "Mass", Mass);
      read(paramtop, "coefs", coefs);
      read(paramtop, "N5", N5);
    }
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception : " << e << endl;
    }
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, EvenOddPrecKNOFermActArrayParams& param)
  {
    EvenOddPrecKNOFermActArrayParams tmp(xml, path);
    param = tmp;
  }



  //! Check stuff
  void EvenOddPrecKNOFermActArray::init()
  {
  }

  void EvenOddPrecKNOFermActArray::initCoeffs(multi1d<Real>& b5_arr,
					      multi1d<Real>& c5_arr) const
  {
    b5_arr.resize(N5);
    c5_arr.resize(N5);


    QDPIO::cout << "Initing General NEF Linop: N5=" << N5 <<endl ;
    QDPIO::cout << "                           a5=" << a5 <<endl ;
    for(int i = 0; i < N5; i++)
      QDPIO::cout<<"                           coef("<<i<<") = "<<coefs[i]<<endl ;

    
    for(int i = 0; i < N5; i++) { 
      b5_arr[i] = Real(0.5)*( coefs[i]  + a5 );
      c5_arr[i] = Real(0.5)*( coefs[i] -  a5 );
    
      QDPIO::cout << " b5["<< i << "] ="<< b5_arr[i] 
		  << " c5["<< i << "] ="<< c5_arr[i] << endl;
    }

  }

  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const EvenOddPrecDWLinOpBaseArray<LatticeFermion>* 
  EvenOddPrecKNOFermActArray::linOp(Handle<const ConnectState> state) const
  {
    multi1d<Real> b5_arr;
    multi1d<Real> c5_arr;

    // Cast the state up to an overlap state
    initCoeffs(b5_arr,c5_arr);
    
    return new EvenOddPrecGenNEFDWLinOpArray(state->getLinks(),OverMass,b5_arr,c5_arr,Mass,N5);
  }

  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator<multi1d<LatticeFermion> >* 
  EvenOddPrecKNOFermActArray::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
  }

  //! Produce a linear operator for this action but with quark mass 1
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const UnprecDWLinOpBaseArray<LatticeFermion>* 
  EvenOddPrecKNOFermActArray::linOpPV(Handle<const ConnectState> state) const
  {
    multi1d<Real> b5_arr;
    multi1d<Real> c5_arr;

    // Cast the state up to an overlap state
    initCoeffs(b5_arr,c5_arr);
    
    return new UnprecNEFDWLinOpArray(state->getLinks(),OverMass,b5_arr,c5_arr,1.0,N5);  // fixed to quark mass 1
  }

  //! Produce an unpreconditioned linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const UnprecDWLinOpBaseArray<LatticeFermion>* 
  EvenOddPrecKNOFermActArray::unprecLinOp(Handle<const ConnectState> state) const
  {
    multi1d<Real> b5_arr;
    multi1d<Real> c5_arr;

    // Cast the state up to an overlap state
    initCoeffs(b5_arr,c5_arr);
    
    return new UnprecNEFDWLinOpArray(state->getLinks(),OverMass,b5_arr,c5_arr,Mass,N5);  // fixed to quark mass 1
  }

  //! Create a ConnectState with just the gauge fields
  const OverlapConnectState*
  EvenOddPrecKNOFermActArray::createState(const multi1d<LatticeColorMatrix>& u_) const
  {
    const OverlapConnectState *ret_val;
    try { 
      ret_val = 
	OverlapConnectStateEnv::createOverlapState(u_, 
						   getFermBC()
						   );
    } 
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
    
    return ret_val;
  }
  
  //! Create a ConnectState with just the gauge fields, and a lower
  //  approximation bound
  const OverlapConnectState*
  EvenOddPrecKNOFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				      const Real& approxMin_) const 
  {
    const OverlapConnectState *ret_val;
    try { 
      ret_val = 
	OverlapConnectStateEnv::createOverlapState(u_, 
						   getFermBC(),
						   approxMin_
						   );
    } 
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }

    return ret_val;
  }

  //! Create a connect State with just approximation range bounds
  const OverlapConnectState*
  EvenOddPrecKNOFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				      const Real& approxMin_,
				      const Real& approxMax_) const
  {
    const OverlapConnectState *ret_val;
    try { 
      ret_val = 
	OverlapConnectStateEnv::createOverlapState(u_, 
						   getFermBC(),
						   approxMin_,
						   approxMax_
						   );
    } 
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }

    return ret_val;
  }
  
  //! Create OverlapConnectState with eigenvalues/vectors
  const OverlapConnectState*
  EvenOddPrecKNOFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				      const multi1d<Real>& lambda_lo_, 
				      const multi1d<LatticeFermion>& evecs_lo_,
				      const Real& lambda_hi_) const
  {
    const OverlapConnectState *ret_val;
    try { 
      ret_val = 
	OverlapConnectStateEnv::createOverlapState(u_, 
						   getFermBC(),
						   lambda_lo_, 
						   evecs_lo_, 
						   lambda_hi_);
      
    } 
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
    
    return ret_val;
  }    
  
  //! Create OverlapConnectState from XML
  const OverlapConnectState*
  EvenOddPrecKNOFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				      XMLReader& state_info_xml,
				      const string& state_info_path) const
  {
    multi1d<LatticeColorMatrix> u_tmp = u_;
    
    // HACK UP A LINEAR OPERATOR TO CHECK EIGENVALUES/VECTORS WITH
    getFermBC().modifyU(u_tmp);
    Handle< const ConnectState > state_aux = new SimpleConnectState(u_tmp);
    Handle< FermBC<LatticeFermion> > fbc_aux = new PeriodicFermBC<LatticeFermion>;
    Real aux_over_mass = -OverMass;
    UnprecWilsonFermAct S_aux( fbc_aux, aux_over_mass );
    
    Handle< const LinearOperator<LatticeFermion> > Maux = 
      S_aux.gamma5HermLinOp(state_aux);
    
    const OverlapConnectState *ret_val;
    
    try {
      ret_val = 
	OverlapConnectStateEnv::createOverlapState(
						   u_,
						   getFermBC(),
						   state_info_xml,
						   state_info_path,
						   *Maux);
    }
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
    
    return ret_val;
  }

  //! Create OverlapConnectState from XML
  const OverlapConnectState*
  EvenOddPrecKNOFermActArray::createState(const multi1d<LatticeColorMatrix>& u_,
				      const OverlapStateInfo& state_info) const
  {

    multi1d<LatticeColorMatrix> u_tmp = u_;
    getFermBC().modifyU(u_tmp);
    Handle< const ConnectState > state_aux = new SimpleConnectState(u_tmp);
    Handle< FermBC<LatticeFermion> > fbc_aux = new PeriodicFermBC<LatticeFermion>;

    Real aux_over_mass = -OverMass;
    UnprecWilsonFermAct S_aux( fbc_aux, aux_over_mass );
    
    Handle< const LinearOperator<LatticeFermion> > Maux = 
      S_aux.gamma5HermLinOp(state_aux);
    
    const OverlapConnectState* ret_val;    
    try {
      ret_val = 
	OverlapConnectStateEnv::createOverlapState(
						   u_,
						   getFermBC(),
						   state_info,
						   *Maux);
    }
    catch(const string& e) { 
      QDPIO::cerr << "Caught Exception: " << e << endl;
      QDP_abort(1);
    }
    
    return ret_val;

  }


}
