// $Id: unprec_parwilson_fermact_w.cc,v 1.5 2004-12-24 04:23:20 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action with parity breaking term
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_parwilson_fermact_w.h"
#include "actions/ferm/linop/unprec_parwilson_linop_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace UnprecParWilsonFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion>* createFermAct(XMLReader& xml_in,
						     const std::string& path)
    {
      return new UnprecParWilsonFermAct(WilsonTypeFermBCEnv::reader(xml_in, path), 
					UnprecParWilsonFermActParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "PARWILSON";

    //! Register the ParWilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct); 
  }


  //! Read parameters
  UnprecParWilsonFermActParams::UnprecParWilsonFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "Mass", Mass);
    read(paramtop, "H", H);
  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecParWilsonFermActParams& param)
  {
    UnprecParWilsonFermActParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >*
  UnprecParWilsonFermAct::linOp(Handle<const ConnectState> state) const
  {
    return new UnprecParWilsonLinOp(state->getLinks(),param.Mass,param.H); 
  }

  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state    gauge field     	       (Read)
   */
  const LinearOperator<LatticeFermion>*
  UnprecParWilsonFermAct::lMdagM(Handle<const ConnectState> state) const
  {
  return new lmdagm<LatticeFermion>(linOp(state));
  }

}
