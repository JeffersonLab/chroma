// $Id: eoprec_dwf_fermact_array_w.cc,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/eoprec_dwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_dwf_linop_array_w.h"
#include "actions/ferm/linop/eoprec_dwf_linop_array_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

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
    WilsonTypeFermAct5D<LatticeFermion, 
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
								       const std::string& path)
    {
      return new EvenOddPrecDWFermActArray(CreateFermStateEnv::reader(xml_in, path), 
					   EvenOddPrecDWFermActArrayParams(xml_in, path));
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
    const std::string name = "DWF";

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


  // Ancient obsolete crap
  namespace
  {

    //! Parameters for chiral fermion actions
    /*! Ancient obsolete junk */
    struct ChiralParam_t
    {
      ChiralParam_t();  // default constructor
      ~ChiralParam_t() {}

      Real       OverMass;
      int        N5;
      Real       a5;
      int        NWilsVec;
    };

    //! Read chiral action like parameters
    void read(XMLReader& xml, const string& path, ChiralParam_t& param);

    //! Write chiral action like parameters
    void write(XMLWriter& xml, const string& path, const ChiralParam_t& param);

    //! Initialize a chiral param struct
    ChiralParam_t::ChiralParam_t()
    {
      OverMass = 0;
      N5       = 0;
      a5       = 1;
      NWilsVec = 0;
    }

    //! Read chiral action like parameters
    void read(XMLReader& xml, const string& path, ChiralParam_t& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "OverMass", param.OverMass);
      read(paramtop, "N5", param.N5);

      string xpath;
      xpath = "a5";
      if (paramtop.count(xpath) != 0)
	read(paramtop, xpath, param.a5);
      else
	param.a5 = 1;

      xpath = "NWilsVec";
      if (paramtop.count(xpath) != 0)
	read(paramtop, xpath, param.NWilsVec);
      else
	param.NWilsVec = 0;
    }

    //! Write chiral action like parameters
    void write(XMLWriter& xml, const string& path, const ChiralParam_t& param)
    {
      push(xml, path);

      write(xml, "OverMass", param.OverMass);
      write(xml, "N5", param.N5);
      write(xml, "a5", param.a5);
      write(xml, "NWilsVec", param.NWilsVec);

      pop(xml);
    }

  } // end namespace



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
  EvenOddPrecDWLikeLinOpBaseArray<LatticeFermion, 
				  multi1d<LatticeColorMatrix>,
				  multi1d<LatticeColorMatrix> >*
  EvenOddPrecDWFermActArray::precLinOp(Handle< FermState<T,P,Q> > state,
				       const Real& m_q) const
  {
    return new EvenOddPrecDWLinOpArray(state,param.OverMass,m_q,param.N5,param.anisoParam);
  }

  //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
  UnprecDWLikeLinOpBaseArray<LatticeFermion, 
			     multi1d<LatticeColorMatrix>,
			     multi1d<LatticeColorMatrix> >*
  EvenOddPrecDWFermActArray::unprecLinOp(Handle< FermState<T,P,Q> > state,
					 const Real& m_q) const
  {
    return new UnprecDWLinOpArray(state,param.OverMass,m_q,param.N5,param.anisoParam);
  }
  

  // Return possibly optimized quark prop solver, solution of preconditioned system
  SystemSolverArray<LatticeFermion>* 
  EvenOddPrecDWFermActArray::qpropT(Handle< FermState<T,P,Q> > state, 
				    const GroupXML_t& invParam) const
  {
    return new DWFQpropT(linOp(state), invLinOp(state,invParam), state, 
			 param.OverMass, param.Mass, param.anisoParam, invParam);
  }


  // Given a complete propagator as a source, this does all the inversions needed
  void 
  EvenOddPrecDWFermActArray::quarkProp(LatticePropagator& q_sol, 
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
      Handle< SystemSolverArray<LatticeFermion> > qpropT(this->qpropT(state, invParam));
      dwf_quarkProp4(q_sol, xml_out, q_src, t_src, j_decay, 
		     qpropT, state, getQuarkMass(), ncg_had);
    }
    else
    {
      Handle< SystemSolver<LatticeFermion> > qprop(this->qprop(state,invParam));
      quarkProp4(q_sol, xml_out, q_src, qprop, quarkSpinType, ncg_had);
    }
  }

}
