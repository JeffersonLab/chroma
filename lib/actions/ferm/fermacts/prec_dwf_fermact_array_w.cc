// $Id: prec_dwf_fermact_array_w.cc,v 1.6 2004-09-08 02:48:25 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_dwf_linop_array_w.h"
#include "actions/ferm/linop/prec_dwf_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"

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

    //! Name to be used
    const std::string name = "DWF";

    //! Register the Wilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActArrayFactory::Instance().registerObject(name, createFermAct); 
  }


  //! Read parameters
  EvenOddPrecDWFermActArrayParams::EvenOddPrecDWFermActArrayParams(XMLReader& xml, 
								   const std::string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    read(paramtop, "WilsonMass", WilsonMass);
    read(paramtop, "m_q", m_q);
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



  //! Produce a linear operator for this action
  /*!
   * \ingroup fermact
   *
   * The operator acts on the odd sublattice
   *
   * \param state 	    gauge field     	       (Read)
   */
  const EvenOddPrecLinearOperator<multi1d<LatticeFermion> >*
  EvenOddPrecDWFermActArray::linOp(Handle<const ConnectState> state) const
  {
    return new EvenOddPrecDWLinOpArray(state->getLinks(),WilsonMass,m_q,N5);
  }

  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the odd sublattice
   *
   * \param state 	    gauge field     	       (Read)
   */
  const LinearOperator<multi1d<LatticeFermion> >*
  EvenOddPrecDWFermActArray::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
  }

  //! Produce a linear operator for this action but with quark mass 1
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator<multi1d<LatticeFermion> >*
  EvenOddPrecDWFermActArray::linOpPV(Handle<const ConnectState> state) const
  {
    // For the PV operator, use the **unpreconditioned** one
    // fixed to quark mass 1
    return new UnprecDWLinOpArray(state->getLinks(),WilsonMass,1.0,N5);
  }

}
