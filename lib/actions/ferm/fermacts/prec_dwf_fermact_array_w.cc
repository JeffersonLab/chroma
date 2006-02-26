// $Id: prec_dwf_fermact_array_w.cc,v 2.2 2006-02-26 03:47:51 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_dwf_linop_array_w.h"
#include "actions/ferm/linop/prec_dwf_linop_array_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

#include "actions/ferm/qprop/quarkprop4_w.h"
#include "actions/ferm/qprop/dwf_quarkprop4_w.h"
#include "io/param_io.h"

#include "actions/ferm/qprop/dwf_qpropt_w.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecDWFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
											const std::string& path)
    {
      return new EvenOddPrecDWFermActArray(WilsonTypeFermBCArrayEnv::reader(xml_in, path), 
					   EvenOddPrecDWFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion>* createFermAct(XMLReader& xml_in,
						 const std::string& path)
    {
      return createFermAct5D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "DWF";

    //! Register all the factories
    bool registerAll()
    {
      return Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct)
	   & Chroma::TheWilsonTypeFermAct5DFactory::Instance().registerObject(name, createFermAct5D);
    }

    //! Register the fermact
    const bool registered = registerAll();
  }


  //! Read parameters
  EvenOddPrecDWFermActArrayParams::EvenOddPrecDWFermActArrayParams(XMLReader& xml, 
								   const std::string& path)
  {
    XMLReader paramtop(xml, path);

    // Check for ancient obsolete tags - try to maintain some backwards compatibility
    if (paramtop.count("ChiralParam") != 0)
    {
      ChiralParam_t pp;
      read(paramtop, "ChiralParam", pp);
      OverMass = pp.OverMass;
      N5 = pp.N5;
      a5 = pp.a5;
    }
    else
    {
      // Read the stuff for the action
      read(paramtop, "OverMass", OverMass);
      read(paramtop, "Mass", Mass);
      read(paramtop, "N5", N5);

      if (paramtop.count("a5") != 0) 
	read(paramtop, "a5", a5);
      else
	a5 = 1.0;
    }

    //  Read optional anisoParam.
    if (paramtop.count("AnisoParam") != 0) 
      read(paramtop, "AnisoParam", anisoParam);
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, EvenOddPrecDWFermActArrayParams& param)
  {
    EvenOddPrecDWFermActArrayParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a preconditioned linear operator for this action with arbitrary quark mass
  const EvenOddPrecDWLikeLinOpBaseArray< LatticeFermion, multi1d<LatticeColorMatrix> >*
  EvenOddPrecDWFermActArray::precLinOp(Handle<const ConnectState> state,
				       const Real& m_q) const
  {
    return new EvenOddPrecDWLinOpArray(state->getLinks(),param.OverMass,m_q,param.N5,param.anisoParam);
  }

  //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
  const UnprecDWLikeLinOpBaseArray< LatticeFermion, multi1d<LatticeColorMatrix> >*
  EvenOddPrecDWFermActArray::unprecLinOp(Handle<const ConnectState> state,
					 const Real& m_q) const
  {
    return new UnprecDWLinOpArray(state->getLinks(),param.OverMass,m_q,param.N5,param.anisoParam);
  }
  

  // Return possibly optimized quark prop solver, solution of preconditioned system
  const SystemSolver< multi1d<LatticeFermion> >* EvenOddPrecDWFermActArray::qpropT(
    Handle<const ConnectState> state,
    const InvertParam_t& invParam) const
  {
    return new DWFQpropT(linOp(state), state, 
			 param.OverMass, param.Mass, param.anisoParam, invParam);
  }


  // Given a complete propagator as a source, this does all the inversions needed
  void 
  EvenOddPrecDWFermActArray::quarkProp(LatticePropagator& q_sol, 
				       XMLWriter& xml_out,
				       const LatticePropagator& q_src,
				       int t_src, int j_decay,
				       Handle<const ConnectState> state,
				       const InvertParam_t& invParam,
				       QuarkSpinType quarkSpinType,
				       bool obsvP,
				       int& ncg_had)
  {
    if (obsvP)
    {
      Handle< const SystemSolver< multi1d<LatticeFermion> > > qpropT(this->qpropT(state,invParam));
      dwf_quarkProp4(q_sol, xml_out, q_src, t_src, j_decay, 
		     qpropT, state, getQuarkMass(), ncg_had);
    }
    else
    {
      Handle< const SystemSolver<LatticeFermion> > qprop(this->qprop(state,invParam));
      quarkProp4(q_sol, xml_out, q_src, qprop, quarkSpinType, ncg_had);
    }
  }

}
