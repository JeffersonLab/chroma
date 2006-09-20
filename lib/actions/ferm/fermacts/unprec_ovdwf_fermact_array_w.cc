// $Id: unprec_ovdwf_fermact_array_w.cc,v 3.2 2006-09-20 20:27:59 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Overlap-DWF (Borici) action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_ovdwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_ovdwf_linop_array_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace UnprecOvDWFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion,
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
								       const std::string& path)
    {
      return new UnprecOvDWFermActArray(CreateFermStateEnv::reader(xml_in, path), 
					UnprecOvDWFermActArrayParams(xml_in, path));
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
    const std::string name = "UNPRECONDITIONED_OVDWF";

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
  UnprecOvDWFermActArrayParams::UnprecOvDWFermActArrayParams(XMLReader& xml, 
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
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecOvDWFermActArrayParams& param)
  {
    UnprecOvDWFermActArrayParams tmp(xml, path);
    param = tmp;
  }

  
  //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
  UnprecDWLikeLinOpBaseArray<LatticeFermion,
			     multi1d<LatticeColorMatrix>,
			     multi1d<LatticeColorMatrix> >* 
  UnprecOvDWFermActArray::unprecLinOp(Handle< FermState<T,P,Q> > state, 
				      const Real& m_q) const
  {
    return new UnprecOvDWLinOpArray(state,OverMass,m_q,N5);
  }

}

