// $Id: prec_nef_fermact_array_w.cc,v 1.8 2004-09-19 02:39:45 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned NEF fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_nef_fermact_array_w.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"
#include "actions/ferm/linop/prec_nef_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/fermfactory_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecNEFFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct< multi1d<LatticeFermion> >* createFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
								XMLReader& xml_in,
								const std::string& path)
    {
      return new EvenOddPrecNEFFermActArray(fbc, EvenOddPrecNEFFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    EvenOddPrecDWFermActBaseArray<LatticeFermion>* createDWFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
								   XMLReader& xml_in,
								   const std::string& path)
    {
      return new EvenOddPrecNEFFermActArray(fbc, EvenOddPrecNEFFermActArrayParams(xml_in, path));
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



  //! Produce a linear operator for this action
  /*!
   * \ingroup fermact
   *
   * The operator acts on the odd sublattice
   *
   * \param state 	    gauge field     	       (Read)
   */
  const EvenOddPrecDWLinOpBaseArray<LatticeFermion>*
  EvenOddPrecNEFFermActArray::linOp(Handle<const ConnectState> state) const
  {
    return new EvenOddPrecNEFDWLinOpArray(state->getLinks(),OverMass,b5,c5,Mass,N5);
  }

  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the odd sublattice
   *
   * \param state 	    gauge field     	       (Read)
   */
  const LinearOperator<multi1d<LatticeFermion> >*
  EvenOddPrecNEFFermActArray::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
  }

  //! Produce a linear operator for this action but with quark mass 1
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const UnprecDWLinOpBaseArray<LatticeFermion>*
  EvenOddPrecNEFFermActArray::linOpPV(Handle<const ConnectState> state) const
  {
    // For the PV operator, use the **unpreconditioned** one
    // fixed to quark mass 1
    return new UnprecNEFDWLinOpArray(state->getLinks(),OverMass,b5,c5,1.0,N5);
  }

}

