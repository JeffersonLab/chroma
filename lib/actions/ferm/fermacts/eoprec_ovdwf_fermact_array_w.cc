// $Id: eoprec_ovdwf_fermact_array_w.cc,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned Overlap-DWF (Borici) action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/eoprec_ovdwf_fermact_array_w.h"
#include "actions/ferm/linop/unprec_ovdwf_linop_array_w.h"
#include "actions/ferm/linop/eoprec_ovdwf_linop_array_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecOvDWFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion, 
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
								       const std::string& path)
    {
      return new EvenOddPrecOvDWFermActArray(CreateFermStateEnv::reader(xml_in, path), 
					     EvenOddPrecOvDWFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* createFermAct(XMLReader& xml_in,
							       const std::string& path)
    {
      return createFermAct5D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "OVDWF"; 

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
	success &= Chroma::TheWilsonTypeFermAct5DFactory::Instance().registerObject(name, createFermAct5D);
	registered = true;
      }
      return success;
    }
  }


  //! Read parameters
  EvenOddPrecOvDWFermActArrayParams::EvenOddPrecOvDWFermActArrayParams(XMLReader& xml, 
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
  void read(XMLReader& xml, const string& path, EvenOddPrecOvDWFermActArrayParams& param)
  {
    EvenOddPrecOvDWFermActArrayParams tmp(xml, path);
    param = tmp;
  }


  //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
  EvenOddPrecDWLikeLinOpBaseArray<LatticeFermion, 
				  multi1d<LatticeColorMatrix>,
				  multi1d<LatticeColorMatrix> >*
  EvenOddPrecOvDWFermActArray::precLinOp(Handle< FermState<T,P,Q> > state,
					 const Real& m_q) const
  {
    return new EvenOddPrecOvDWLinOpArray(state,OverMass,m_q,N5);
  }

  //! Produce an even-odd preconditioned linear operator for this action with arbitrary quark mass
  UnprecDWLikeLinOpBaseArray<LatticeFermion, 
			     multi1d<LatticeColorMatrix>,
			     multi1d<LatticeColorMatrix> >*
  EvenOddPrecOvDWFermActArray::unprecLinOp(Handle< FermState<T,P,Q> > state,
					   const Real& m_q) const
  {
    return new UnprecOvDWLinOpArray(state,OverMass,m_q,N5);
  }

}

