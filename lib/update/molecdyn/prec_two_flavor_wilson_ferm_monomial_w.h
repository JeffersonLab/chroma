#ifndef PREC_TWO_FLAVOR_WILSON_MONOMIAL_H
#define PREC_TWO_FLAVOR_WILSON_MONOMIAL_H

#include "chromabase.h"

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/abs_monomial.h"
#include "actions/ferm/fermacts/prec_wilson_fermact_w.h"

using namespace std;
using namespace QDP;
using namespace std;

namespace Chroma {

  namespace EvenOddPrecTwoFlavorWilsonFermMonomialEnv {
    extern const string name;
    extern const bool registered;
  };

  // Parameter structure
  struct EvenOddPrecTwoFlavorWilsonFermMonomialParams {
    // Base Constructor
    EvenOddPrecTwoFlavorWilsonFermMonomialParams();

    // Read monomial from some root path
    EvenOddPrecTwoFlavorWilsonFermMonomialParams(XMLReader& in, const std::string&  path);
    InvertParam_t inv_param; // Inverter Parameters
    string ferm_act;
  };

  void read(XMLReader& xml, const string& path, EvenOddPrecTwoFlavorWilsonFermMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const EvenOddPrecTwoFlavorWilsonFermMonomialParams& params);


  class EvenOddPrecTwoFlavorWilsonFermMonomial :
    public  TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a FermBC and a parameter struct
      EvenOddPrecTwoFlavorWilsonFermMonomial(Handle< FermBC<LatticeFermion> > fbc_,
					const EvenOddPrecTwoFlavorWilsonFermMonomialParams& param_);


      // Construct from a fermact handle and inv params
      // FermAct already holds BC-s
      EvenOddPrecTwoFlavorWilsonFermMonomial(Handle< const EvenOddPrecWilsonFermAct >& fermact_, const InvertParam_t& inv_param_ ) : fermact(fermact_), inv_param(inv_param_) {}

      // Copy Constructor
      EvenOddPrecTwoFlavorWilsonFermMonomial(const EvenOddPrecTwoFlavorWilsonFermMonomial& m) : fermact((m.fermact)), inv_param(m.inv_param) 
	{
	  phi = m.phi;
	}


      const LatticeFermion& debugGetPhi(void) const {
	return getPhi();
      }

      void debugGetX(LatticeFermion& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const {
	getX(X,s);
      }

      const EvenOddPrecWilsonFermAct& debugGetFermAct(void) const { 
	return getFermAct();
      }
      
 
      //! Even even contribution (eg ln det Clover)
      Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
			                     multi1d<LatticeColorMatrix> >& s) const {
	return Double(0);
      }


    protected:


      LatticeFermion& getPhi(void) {
	return phi;
      }

      const LatticeFermion& getPhi(void) const {
	return phi;
      }

      const EvenOddPrecWilsonFermAct& getFermAct(void) const { 
	return *fermact;
      }

      // Do inversion M^dag M X = phi
      void getX(LatticeFermion& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const;


      
    private:
 
      // Pseudofermion field phi
      LatticeFermion phi;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const EvenOddPrecWilsonFermAct> fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;
    };


}; //end namespace chroma













#endif
