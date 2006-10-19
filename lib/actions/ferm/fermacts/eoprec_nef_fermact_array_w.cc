// $Id: eoprec_nef_fermact_array_w.cc,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned NEF fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/eoprec_nef_fermact_array_w.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"
#include "actions/ferm/linop/eoprec_nef_linop_array_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "actions/ferm/qprop/quarkprop4_w.h"
#include "actions/ferm/qprop/nef_quarkprop4_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecNEFFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion, 
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
								       const std::string& path)
    {
      return new EvenOddPrecNEFFermActArray(CreateFermStateEnv::reader(xml_in, path), 
					    EvenOddPrecNEFFermActArrayParams(xml_in, path));
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
    const std::string name = "NEF";

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
  EvenOddPrecNEFFermActArrayParams::EvenOddPrecNEFFermActArrayParams(XMLReader& xml, 
							     const std::string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    read(paramtop, "OverMass", OverMass);
    read(paramtop, "Mass", Mass);
    read(paramtop, "N5", N5);
    read(paramtop, "b5", b5);
    read(paramtop, "c5", c5);
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, EvenOddPrecNEFFermActArrayParams& param)
  {
    EvenOddPrecNEFFermActArrayParams tmp(xml, path);
    param = tmp;
  }



  //! Produce an even-odd preconditioned linear operator for this action with arbitrary quark mass
  EvenOddPrecDWLikeLinOpBaseArray<LatticeFermion, 
				  multi1d<LatticeColorMatrix>,
				  multi1d<LatticeColorMatrix> >*
  EvenOddPrecNEFFermActArray::precLinOp(Handle< FermState<T,P,Q> > state,
					const Real& m_q) const
  {
    return new EvenOddPrecNEFDWLinOpArray(state,params.OverMass,
					  params.b5,params.c5,m_q,params.N5);
  }

  //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
  UnprecDWLikeLinOpBaseArray<LatticeFermion, 
			     multi1d<LatticeColorMatrix>,
			     multi1d<LatticeColorMatrix> >*
  EvenOddPrecNEFFermActArray::unprecLinOp(Handle< FermState<T,P,Q> > state,
					  const Real& m_q) const
  {
    multi1d<Real> bb5(params.N5);
    multi1d<Real> cc5(params.N5);

    bb5 = params.b5;
    cc5 = params.c5;

    return new UnprecNEFDWLinOpArray(state,params.OverMass,bb5,cc5,m_q,params.N5);
  }

  
  // Given a complete propagator as a source, this does all the inversions needed
  void 
  EvenOddPrecNEFFermActArray::quarkProp(LatticePropagator& q_sol, 
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

