// $Id: unprec_nef_fermact_array_w.cc,v 1.9 2004-10-03 01:21:19 edwards Exp $
/*! \file
 *  \brief Unpreconditioned NEF fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_nef_fermact_array_w.h"
#include "actions/ferm/linop/unprec_nef_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/fermfactory_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace UnprecNEFFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct< multi1d<LatticeFermion> >* createFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
								XMLReader& xml_in,
								const std::string& path)
    {
      return new UnprecNEFFermActArray(fbc, UnprecNEFFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    UnprecDWFermActBaseArray<LatticeFermion>* createDWFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
							      XMLReader& xml_in,
							      const std::string& path)
    {
      return new UnprecNEFFermActArray(fbc, UnprecNEFFermActArrayParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_NEF";

    //! Register the Wilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActArrayFactory::Instance().registerObject(name, createFermAct)
                          & Chroma::TheUnprecDWFermActBaseArrayFactory::Instance().registerObject(name, createDWFermAct); 
  }


  //! Read parameters
  UnprecNEFFermActArrayParams::UnprecNEFFermActArrayParams(XMLReader& xml, 
							   const std::string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    read(paramtop, "OverMass", OverMass);
    read(paramtop, "Mass", Mass);
    read(paramtop, "N5", N5);

    b5.resize(N5);
    c5.resize(N5);
    multi1d<Real> bb5;
    multi1d<Real> cc5;

    read(paramtop, "b5", bb5);
    read(paramtop, "c5", cc5);

    if (bb5.size() != cc5.size())
    {
      QDPIO::cerr << "UnprecNEF: array sizes incorrect" << endl;
      QDP_abort(1);
    }

    if (bb5.size() != N5)
    {
      if (bb5.size() == 1)
      {
	for(int n=0; n < N5; ++n)
	{
	  b5[n] = bb5[0];
	  c5[n] = cc5[0];
	}
      }
      else
      {
	QDPIO::cerr << "UnprecNEF: array sizes incorrect" << endl;
	QDP_abort(1);
      }
    }
    else
    {
      b5 = bb5;
      c5 = cc5;
    }
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecNEFFermActArrayParams& param)
  {
    UnprecNEFFermActArrayParams tmp(xml, path);
    param = tmp;
  }



  //! Check stuff
  void UnprecNEFFermActArray::init()
  {
    if ((b5.size() != c5.size()) || (b5.size() != N5))
    {
      QDPIO::cerr << "UnprecNEF: array sizes incompatible" << endl;
      QDP_abort(1);
    }
  }

  //! Produce a linear operator for this action
  /*!
   * \ingroup fermact
   *
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const UnprecDWLinOpBaseArray<LatticeFermion>* 
  UnprecNEFFermActArray::linOp(Handle<const ConnectState> state) const
  {
    return new UnprecNEFDWLinOpArray(state->getLinks(),OverMass,b5,c5,Mass,N5);
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
  UnprecNEFFermActArray::lMdagM(Handle<const ConnectState> state) const
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
  const UnprecDWLinOpBaseArray<LatticeFermion>* 
  UnprecNEFFermActArray::linOpPV(Handle<const ConnectState> state) const
  {
    return new UnprecNEFDWLinOpArray(state->getLinks(),OverMass,b5,c5,1.0,N5);  // fixed to quark mass 1
  }

}
