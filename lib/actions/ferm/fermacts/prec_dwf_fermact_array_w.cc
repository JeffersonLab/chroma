// $Id: prec_dwf_fermact_array_w.cc,v 1.11 2004-12-09 03:58:03 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_dwf_linop_array_w.h"
#include "actions/ferm/linop/prec_dwf_linop_array_w.h"

#include "actions/ferm/fermacts/fermfactory_w.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecDWFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct< multi1d<LatticeFermion> >* createFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
								XMLReader& xml_in,
								const std::string& path)
    {
      return new EvenOddPrecDWFermActArray(fbc, EvenOddPrecDWFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    EvenOddPrecDWFermActBaseArray<LatticeFermion>* createDWFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
								   XMLReader& xml_in,
								   const std::string& path)
    {
      return new EvenOddPrecDWFermActArray(fbc, EvenOddPrecDWFermActArrayParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "DWF";

    //! Register the Wilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActArrayFactory::Instance().registerObject(name, createFermAct)
                          & Chroma::TheEvenOddPrecDWFermActBaseArrayFactory::Instance().registerObject(name, createDWFermAct); 
  }


  //! Read parameters
  EvenOddPrecDWFermActArrayParams::EvenOddPrecDWFermActArrayParams(XMLReader& xml, 
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
  void read(XMLReader& xml, const string& path, EvenOddPrecDWFermActArrayParams& param)
  {
    EvenOddPrecDWFermActArrayParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a preconditioned linear operator for this action with arbitrary quark mass
  const EvenOddPrecDWLinOpBaseArray<LatticeFermion>*
  EvenOddPrecDWFermActArray::precLinOp(Handle<const ConnectState> state,
				       const Real& m_q) const
  {
    return new EvenOddPrecDWLinOpArray(state->getLinks(),OverMass,m_q,N5);
  }

  //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
  const UnprecDWLinOpBaseArray<LatticeFermion>*
  EvenOddPrecDWFermActArray::unprecLinOp(Handle<const ConnectState> state,
					 const Real& m_q) const
  {
    return new UnprecDWLinOpArray(state->getLinks(),OverMass,m_q,N5);
  }

}
