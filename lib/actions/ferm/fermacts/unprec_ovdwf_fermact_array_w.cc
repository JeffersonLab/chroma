// $Id: unprec_ovdwf_fermact_array_w.cc,v 1.11 2004-12-24 04:23:20 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Overlap-DWF (Borici) action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_ovdwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_ovdwf_linop_array_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace UnprecOvDWFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct< multi1d<LatticeFermion> >* createFermAct(XMLReader& xml_in,
								const std::string& path)
    {
      return new UnprecOvDWFermActArray(WilsonTypeFermBCArrayEnv::reader(xml_in, path), 
					UnprecOvDWFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    UnprecDWFermActBaseArray<LatticeFermion>* createDWFermAct(XMLReader& xml_in,
							      const std::string& path)
    {
      return new UnprecOvDWFermActArray(WilsonTypeFermBCArrayEnv::reader(xml_in, path), 
					UnprecOvDWFermActArrayParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_OVDWF";

    //! Register the Wilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActArrayFactory::Instance().registerObject(name, createFermAct)
                          & Chroma::TheUnprecDWFermActBaseArrayFactory::Instance().registerObject(name, createDWFermAct); 
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
  const UnprecDWLinOpBaseArray<LatticeFermion>* 
  UnprecOvDWFermActArray::unprecLinOp(Handle<const ConnectState> state, 
				      const Real& m_q) const
  {
    return new UnprecOvDWLinOpArray(state->getLinks(),OverMass,m_q,N5);
  }

}

