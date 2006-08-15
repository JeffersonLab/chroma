#include "chromabase.h"
#include "actions/ferm/fermacts/stout_fermstate_w.h"
#include "actions/ferm/fermbcs/fermbcs.h"
#include "update/molecdyn/monomial/fixed_random_ferm_monomial.h"
#include "update/molecdyn/monomial/monomial_factory.h"

namespace Chroma { 

  namespace FixedRandomFermMonomial4DEnv { 

    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) 
    {
      QDPIO::cout << "Create Monomial: " << name << endl;
      
      return new  FixedRandomFermMonomial4D(FixedRandomFermMonomialParams(xml,path));
      
    }
    
    const std::string name("FIXED_RANDOM_X_FERM_MONOMIAL");

    //! Register all the objects
    bool registerAll()
    {
      bool foo = true;
      
      foo &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
      
      return foo;
    }

    const bool registered = registerAll();
  };

  FixedRandomFermMonomial4D::FixedRandomFermMonomial4D(const FixedRandomFermMonomialParams& p){ 
    
    X.resize(Nd);
    for(int mu=0; mu < Nd; mu++) { 
      gaussian(X[mu]);
      reunit(X[mu]);
    }
    
    std::istringstream is( p.fermstate.xml );
    QDPIO::cout << "FermStateXML is: " << p.fermstate.xml << endl;
    
    XMLReader xml_in(is);
    
    cfs =  new CreateStoutFermState(WilsonTypeFermBCEnv::reader(xml_in, "/FermState"),
				    StoutFermStateParams(xml_in, "/FermState"));


  }

  void FixedRandomFermMonomial4D::dsdq(P& F, const AbsFieldState<P,Q>& s)
  {
    Handle< StoutFermState > state = (*cfs)(s.getQ()); // Create a 
    F.resize(Nd);
    F = X;
    
    state->deriv(F);
     
  }

  Double  FixedRandomFermMonomial4D::S(const AbsFieldState<P,Q>& s)
  {
    Handle< StoutFermState > state = (*cfs)(s.getQ()); // Create a 
    const Q& u = state->getLinks();
    
    Double ret_val = sum(real(trace(u[0]*X[0])));
    for(int mu=1; mu < Nd; mu++) { 
      ret_val += sum(real(trace(u[mu]*X[mu])));
    }
    
    return Double(2)*ret_val;
  }

};

