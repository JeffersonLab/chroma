// -*- C++ -*-
// $Id: fixed_random_ferm_monomial.cc,v 3.9 2008-05-23 18:40:41 edwards Exp $

/*! @file
 * @brief Fixed random monomial
 */

#include "chromabase.h"
#include "actions/ferm/fermstates/stout_fermstate_w.h"
#include "actions/ferm/fermbcs/fermbcs.h"
#include "update/molecdyn/monomial/fixed_random_ferm_monomial.h"
#include "update/molecdyn/monomial/monomial_factory.h"

namespace Chroma { 

  namespace FixedRandomFermMonomial4DEnv 
  { 
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
    {
      return new  FixedRandomFermMonomial4D(FixedRandomFermMonomialParams(xml,path));
      
    }
    
    const std::string name("FIXED_RANDOM_X_FERM_MONOMIAL");

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
	registered = true;
      }
      return success;
    }
  }


  FixedRandomFermMonomial4D::FixedRandomFermMonomial4D(const FixedRandomFermMonomialParams& p)
  { 
    START_CODE();
   
    X.resize(Nd);
    for(int mu=0; mu < Nd; mu++) { 
      gaussian(X[mu]);
      reunit(X[mu]);
    }
    
    std::istringstream is( p.fermstate.xml );
    QDPIO::cout << "FermStateXML is: " << p.fermstate.xml << endl;
    
    XMLReader xml_in(is);
    
    cfs =  new CreateStoutFermState<LatticeFermion,
      multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >(WilsonTypeFermBCEnv::reader(xml_in, "/FermState"),
				    StoutFermStateParams(xml_in, "/FermState"));

    END_CODE();
  }

  void FixedRandomFermMonomial4D::dsdq(P& F, const AbsFieldState<P,Q>& s)
  {
    START_CODE();

    Handle< StoutFermState<LatticeFermion,
      multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> > 
      > state = (*cfs)(s.getQ()); // Create a
 
    F.resize(Nd);
    F = X;
    
    state->deriv(F);

    XMLWriter& xml_out = TheXMLLogWriter::Instance();

    push(xml_out, "FixedRandomFermMonomial4D"); 
    Double F_sq = norm2(F);
    write(xml_out, "F_sq", F_sq);
    pop(xml_out);
    
    END_CODE();
  }

  Double  FixedRandomFermMonomial4D::S(const AbsFieldState<P,Q>& s)
  {
    START_CODE();

    XMLWriter& xml_out = TheXMLLogWriter::Instance();
    push(xml_out, "FixedRandomFermMonomial4D");

    Handle< StoutFermState<LatticeFermion,
      multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >
      > state = (*cfs)(s.getQ()); // Create a 
    const Q& u = state->getLinks();
    
    Double ret_val = sum(real(trace(u[0]*X[0])));
    for(int mu=1; mu < Nd; mu++) { 
      ret_val += sum(real(trace(u[mu]*X[mu])));
    }
    ret_val *= Double(2);
    write(xml_out, "S", ret_val);
    pop(xml_out);

    return ret_val;
    
    END_CODE();
  }

}

