#ifndef UNPREC_TWO_FLAVOR_WILSON_TYPE_FERM_MONOMIAL_H
#define UNPREC_TWO_FLAVOR_WILSON_TYOE_FERM_MONOMIAL_H

#include "chromabase.h"

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/abs_monomial.h"
//#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "fermact.h"

using namespace std;
using namespace QDP;
using namespace std;

namespace Chroma {

  namespace UnprecTwoFlavorWilsonTypeFermMonomialEnv {
    //extern const string name;
    extern const bool registered;
  };

  // Parameter structure
  struct UnprecTwoFlavorWilsonTypeFermMonomialParams {
    // Base Constructor
    UnprecTwoFlavorWilsonTypeFermMonomialParams();

    // Read monomial from some root path
    UnprecTwoFlavorWilsonTypeFermMonomialParams(XMLReader& in, const std::string&  path);
    InvertParam_t inv_param; // Inverter Parameters
    string ferm_act;
  };

  void read(XMLReader& xml, const string& path, UnprecTwoFlavorWilsonTypeFermMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const UnprecTwoFlavorWilsonTypeFermMonomialParams& params);


  class UnprecTwoFlavorWilsonTypeFermMonomial :
    public  TwoFlavorExactUnprecWilsonTypeFermMonomial< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a FermBC and a parameter struct
      UnprecTwoFlavorWilsonTypeFermMonomial(const UnprecTwoFlavorWilsonTypeFermMonomialParams& param_);


      // Construct from a fermact handle and inv params
      // FermAct already holds BC-s
      UnprecTwoFlavorWilsonTypeFermMonomial(Handle< const UnprecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix> >  >& fermact_, const InvertParam_t& inv_param_ ) : fermact(fermact_), inv_param(inv_param_) {}

      // Copy Constructor
      UnprecTwoFlavorWilsonTypeFermMonomial(const UnprecTwoFlavorWilsonTypeFermMonomial& m) : fermact((m.fermact)), inv_param(m.inv_param) 
	{
	  phi = m.phi;
	}


      const LatticeFermion& debugGetPhi(void) const {
	return getPhi();
      }

      void debugGetX(LatticeFermion& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const {
	getX(X,s);
      }

      const UnprecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix> >& debugGetFermAct(void) const { 
	return getFermAct();
      }
      
 

    protected:


      LatticeFermion& getPhi(void) {
	return phi;
      }

      const LatticeFermion& getPhi(void) const {
	return phi;
      }

      const UnprecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
      }

      // Do inversion M^dag M X = phi
      void getX(LatticeFermion& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const;


      
    private:
 
      // Pseudofermion field phi
      LatticeFermion phi;

      // A handle for the UnprecWilsonTypeFermAct
      Handle<const UnprecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;
    };


}; //end namespace chroma













#endif
