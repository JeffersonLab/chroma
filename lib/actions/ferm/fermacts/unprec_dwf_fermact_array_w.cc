// $Id: unprec_dwf_fermact_array_w.cc,v 1.11 2004-12-09 03:58:03 edwards Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_dwf_linop_array_w.h"

#include "actions/ferm/fermacts/fermfactory_w.h"

using namespace Chroma;

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace UnprecDWFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct< multi1d<LatticeFermion> >* createFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
								XMLReader& xml_in,
								const std::string& path)
    {
      return new UnprecDWFermActArray(fbc, UnprecDWFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Only differs in return type */
    UnprecDWFermActBaseArray<LatticeFermion>* createDWFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
							      XMLReader& xml_in,
							      const std::string& path)
    {
      return new UnprecDWFermActArray(fbc, UnprecDWFermActArrayParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_DWF";

    //! Register the Wilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActArrayFactory::Instance().registerObject(name, createFermAct)
                          & Chroma::TheUnprecDWFermActBaseArrayFactory::Instance().registerObject(name, createDWFermAct); 
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
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecDWFermActArrayParams& param)
  {
    UnprecDWFermActArrayParams tmp(xml, path);
    param = tmp;
  }


  
  //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
  const UnprecDWLinOpBaseArray<LatticeFermion>* 
  UnprecDWFermActArray::unprecLinOp(Handle<const ConnectState> state, 
				    const Real& m_q) const
  {
    return new UnprecDWLinOpArray(state->getLinks(),OverMass,m_q,N5);
  }

}
