// $Id: prec_parwilson_fermact_w.cc,v 1.4 2004-12-12 21:22:15 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action with parity breaking term
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_parwilson_fermact_w.h"
#include "actions/ferm/linop/prec_parwilson_linop_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/fermfactory_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecParWilsonFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion>* createFermAct(Handle< FermBC<LatticeFermion> > fbc,
						     XMLReader& xml_in,
						     const std::string& path)
    {
      return new EvenOddPrecParWilsonFermAct(fbc, EvenOddPrecParWilsonFermActParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_PARWILSON";

    //! Register the ParWilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct); 
  }


  //! Read parameters
  EvenOddPrecParWilsonFermActParams::EvenOddPrecParWilsonFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "Mass", Mass);
    read(paramtop, "H", H);
  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, EvenOddPrecParWilsonFermActParams& param)
  {
    EvenOddPrecParWilsonFermActParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a linear operator for this action
  /*!
   * The operator acts on the odd subset
   *
   * \param state 	    gauge field     	       (Read)
   */
  const EvenOddPrecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* 
  EvenOddPrecParWilsonFermAct::linOp(Handle<const ConnectState> state) const
  {
    return new EvenOddPrecParWilsonLinOp(state->getLinks(),param.Mass,param.H);
  }

  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the odd subset
   *
   * \param state 	    gauge field     	       (Read)
   */
  const LinearOperator<LatticeFermion>* 
  EvenOddPrecParWilsonFermAct::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<LatticeFermion>(linOp(state));
  }

}
