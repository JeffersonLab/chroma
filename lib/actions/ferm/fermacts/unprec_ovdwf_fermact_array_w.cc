// $Id: unprec_ovdwf_fermact_array_w.cc,v 1.8 2004-09-18 17:14:45 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Overlap-DWF (Borici) action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_ovdwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_ovdwf_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/fermfactory_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace UnprecOvDWFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct< multi1d<LatticeFermion> >* createFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
								XMLReader& xml_in,
								const std::string& path)
    {
      return new UnprecOvDWFermActArray(fbc, UnprecOvDWFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    UnprecDWFermActBaseArray<LatticeFermion>* createDWFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
							      XMLReader& xml_in,
							      const std::string& path)
    {
      return new UnprecOvDWFermActArray(fbc, UnprecOvDWFermActArrayParams(xml_in, path));
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



  //! Produce a linear operator for this action
  /*!
   * \ingroup fermact
   *
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator<multi1d<LatticeFermion> >* 
  UnprecOvDWFermActArray::linOp(Handle<const ConnectState> state) const
  {
    return new UnprecOvDWLinOpArray(state->getLinks(),OverMass,Mass,N5);
  }


  //! Produce a M^dag.M linear operator for this action
  /*!
   * \ingroup fermact
   *
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator<multi1d<LatticeFermion> >* 
  UnprecOvDWFermActArray::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
  }


  //! Produce a linear operator for this action but with quark mass 1
  /*!
   * \ingroup fermact
   *
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator<multi1d<LatticeFermion> >* 
  UnprecOvDWFermActArray::linOpPV(Handle<const ConnectState> state) const
  {
    return new UnprecOvDWLinOpArray(state->getLinks(),OverMass,1.0,N5);  // fixed to quark mass 1
  }

}

