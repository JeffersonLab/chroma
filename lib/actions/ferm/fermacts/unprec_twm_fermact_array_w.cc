// $Id: unprec_twm_fermact_array_w.cc,v 1.1 2008-11-04 18:42:58 edwards Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall fermion action
 */

#include "chromabase.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"

#include "actions/ferm/fermacts/unprec_dwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_dwf_linop_array_w.h"

// #include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "actions/ferm/qprop/quarkprop4_w.h"
#include "actions/ferm/qprop/dwf_quarkprop4_w.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace UnprecDWFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion,
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
								       const std::string& path)
    {
      return new UnprecDWFermActArray(CreateFermStateEnv::reader(xml_in, path), 
				      UnprecDWFermActArrayParams(xml_in, path));
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
    const std::string name = "UNPRECONDITIONED_DWF";

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
  UnprecDWFermActArrayParams::UnprecDWFermActArrayParams(XMLReader& xml, 
								   const std::string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    read(paramtop, "OverMass", OverMass);
    read(paramtop, "Mass", Mass);
    read(paramtop, "N5", N5);

    if (paramtop.count("a5") != 0) 
      read(paramtop, "a5", a5);
    else
      a5 = 1.0;

    //  Read optional anisoParam.
    if (paramtop.count("AnisoParam") != 0) 
      read(paramtop, "AnisoParam", anisoParam);
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecDWFermActArrayParams& param)
  {
    UnprecDWFermActArrayParams tmp(xml, path);
    param = tmp;
  }


  
  //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
  UnprecDWLikeLinOpBaseArray<LatticeFermion,
			     multi1d<LatticeColorMatrix>,
			     multi1d<LatticeColorMatrix> >* 
  UnprecDWFermActArray::unprecLinOp(Handle< FermState<T,P,Q> > state, 
				    const Real& m_q) const
  {
    return new UnprecDWLinOpArray(state,param.OverMass,m_q,param.N5,param.anisoParam);
  }

  
  // Given a complete propagator as a source, this does all the inversions needed
  void 
  UnprecDWFermActArray::quarkProp(LatticePropagator& q_sol, 
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
    {
      Handle< SystemSolverArray<T> > qpropT(this->qpropT(state,invParam));
      dwf_quarkProp4(q_sol, xml_out, q_src, t_src, j_decay, qpropT, state, getQuarkMass(), ncg_had);
    }
    else
    {
      Handle< SystemSolver<T> > qprop(this->qprop(state,invParam));
      quarkProp4(q_sol, xml_out, q_src, qprop, quarkSpinType, ncg_had);
    }
  }

}
