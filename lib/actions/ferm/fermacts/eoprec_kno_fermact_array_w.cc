// $Id: eoprec_kno_fermact_array_w.cc,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief preconditioned KNO fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/eoprec_kno_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"
#include "actions/ferm/linop/eoprec_nef_general_linop_array_w.h"

#include "actions/ferm/fermacts/zolotarev.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "actions/ferm/qprop/quarkprop4_w.h"
#include "actions/ferm/qprop/nef_quarkprop4_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecKNOFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion, 
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
								       const std::string& path)
    {
      return new EvenOddPrecKNOFermActArray(CreateFermStateEnv::reader(xml_in, path), 
					    EvenOddPrecKNOFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion, 
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* createFermAct(XMLReader& xml_in,
							       const std::string& path)
    {
      return createFermAct5D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "KNO";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
	success &= Chroma::TheWilsonTypeFermAct5DFactory::Instance().registerObject(name, createFermAct5D);
	registered = true;
      }
      return success;
    }
  }


  //! Read parameters
  EvenOddPrecKNOFermActArrayParams::EvenOddPrecKNOFermActArrayParams(XMLReader& xml, 
								     const std::string& path)
  {
    XMLReader paramtop(xml, path);
    try 
    {
      // Read the stuff for the action
      read(paramtop, "OverMass", OverMass);
      read(paramtop, "Mass", Mass);
      read(paramtop, "a5", a5);
      read(paramtop, "coefs", coefs);
      //read(paramtop, "N5", N5);
      N5 = coefs.size();
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


  //! Produce a preconditioned linear operator for this action with arbitrary quark mass
  EvenOddPrecDWLikeLinOpBaseArray<LatticeFermion, 
				  multi1d<LatticeColorMatrix>,
				  multi1d<LatticeColorMatrix> >* 
  EvenOddPrecKNOFermActArray::precLinOp(Handle< FermState<T,P,Q> > state,
					const Real& m_q) const
  {
    multi1d<Real> b5_arr;
    multi1d<Real> c5_arr;

    // Cast the state up to an overlap state
    initCoeffs(b5_arr,c5_arr);
    
    return new EvenOddPrecGenNEFDWLinOpArray(state,OverMass,b5_arr,c5_arr,m_q,N5);
  }

  //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
  UnprecDWLikeLinOpBaseArray<LatticeFermion, 
			     multi1d<LatticeColorMatrix>,
			     multi1d<LatticeColorMatrix> >* 
  EvenOddPrecKNOFermActArray::unprecLinOp(Handle< FermState<T,P,Q> > state,
					  const Real& m_q) const
  {
    multi1d<Real> b5_arr;
    multi1d<Real> c5_arr;

    // Cast the state up to an overlap state
    initCoeffs(b5_arr,c5_arr);
    
    return new UnprecNEFDWLinOpArray(state,OverMass,b5_arr,c5_arr,m_q,N5);
  }

  
  // Given a complete propagator as a source, this does all the inversions needed
  void 
  EvenOddPrecKNOFermActArray::quarkProp(LatticePropagator& q_sol, 
					XMLWriter& xml_out,
					const LatticePropagator& q_src,
					int t_src, int j_decay,
					Handle< FermState<T,P,Q> > state,
					const GroupXML_t& invParam,
					QuarkSpinType quarkSpinType,
					bool obsvP,
					int& ncg_had) const
  {
    if (obsvP && (quarkSpinType == QUARK_SPIN_TYPE_FULL))
      nef_quarkProp4(q_sol, xml_out, q_src, t_src, j_decay, *this, state, invParam, ncg_had);
    else
    {
      Handle< SystemSolver<LatticeFermion> > qprop(this->qprop(state,invParam));
      quarkProp4(q_sol, xml_out, q_src, qprop, quarkSpinType, ncg_had);
    }
  }


}
