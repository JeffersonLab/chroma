// $Id: prec_nef_fermact_array_w.cc,v 1.12 2004-12-24 04:23:20 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned NEF fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_nef_fermact_array_w.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"
#include "actions/ferm/linop/prec_nef_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecNEFFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct< multi1d<LatticeFermion> >* createFermAct(XMLReader& xml_in,
								const std::string& path)
    {
      return new EvenOddPrecNEFFermActArray(WilsonTypeFermBCArrayEnv::reader(xml_in, path), 
					    EvenOddPrecNEFFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    EvenOddPrecDWFermActBaseArray<LatticeFermion>* createDWFermAct(XMLReader& xml_in,
								   const std::string& path)
    {
      return new EvenOddPrecNEFFermActArray(WilsonTypeFermBCArrayEnv::reader(xml_in, path), 
					    EvenOddPrecNEFFermActArrayParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "NEF";

    //! Register the Wilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActArrayFactory::Instance().registerObject(name, createFermAct)
                          & Chroma::TheEvenOddPrecDWFermActBaseArrayFactory::Instance().registerObject(name, createDWFermAct); 
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
  const EvenOddPrecDWLinOpBaseArray<LatticeFermion>*
  EvenOddPrecNEFFermActArray::precLinOp(Handle<const ConnectState> state,
					const Real& m_q) const
  {
    return new EvenOddPrecNEFDWLinOpArray(state->getLinks(),OverMass,b5,c5,m_q,N5);
  }

  //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
  const UnprecDWLinOpBaseArray<LatticeFermion>*
  EvenOddPrecNEFFermActArray::unprecLinOp(Handle<const ConnectState> state,
					  const Real& m_q) const
  {
    multi1d<Real> bb5(N5);
    multi1d<Real> cc5(N5);

    bb5 = b5;
    cc5 = c5;

    return new UnprecNEFDWLinOpArray(state->getLinks(),OverMass,bb5,cc5,m_q,N5);
  }

}

