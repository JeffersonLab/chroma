#ifndef PREC_TWO_FLAVOR_WILSON_MONOMIAL_H
#define PREC_TWO_FLAVOR_WILSON_MONOMIAL_H

#include "chromabase.h"

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/abs_monomial.h"
#include "fermact.h"

using namespace std;
using namespace QDP;
using namespace std;

namespace Chroma {

  namespace EvenOddPrecTwoFlavorWilsonTypeFermMonomialEnv {
    //    extern const string name;
    extern const bool registered;
  };

  // Parameter structure
  struct EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams {
    // Base Constructor
    EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams();

    // Read monomial from some root path
    EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams(XMLReader& in, const std::string&  path);
    InvertParam_t inv_param; // Inverter Parameters
    string ferm_act;
  };

  void read(XMLReader& xml, const string& path, EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams& params);


  class EvenOddPrecTwoFlavorWilsonTypeFermMonomial :
    public  TwoFlavorExactEvenOddPrecWilsonTypeFermMonomial< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
    {
    public: 
      // Construct out of a FermBC and a parameter struct
      EvenOddPrecTwoFlavorWilsonTypeFermMonomial(const EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams& param_);


      // Construct from a fermact handle and inv params
      // FermAct already holds BC-s
      EvenOddPrecTwoFlavorWilsonTypeFermMonomial(
	Handle< 
        const EvenOddPrecWilsonTypeFermAct< LatticeFermion, 
                                            multi1d<LatticeColorMatrix> 
                                          > 
              > fermact_, 
              const InvertParam_t& inv_param_ ) : fermact(fermact_), inv_param(inv_param_) {}

      // Copy Constructor
      EvenOddPrecTwoFlavorWilsonTypeFermMonomial(const EvenOddPrecTwoFlavorWilsonTypeFermMonomial& m) : fermact((m.fermact)), inv_param(m.inv_param) 
	{
	  phi = m.phi;
	}


      const LatticeFermion& debugGetPhi(void) const {
	return getPhi();
      }

      void debugGetX(LatticeFermion& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const {
	getX(X,s);
      }

      const EvenOddPrecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix> >& debugGetFermAct(void) const { 
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

      const EvenOddPrecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix> >& getFermAct(void) const { 
	return *fermact;
      }

      // Do inversion M^dag M X = phi
      void getX(LatticeFermion& X, const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const;


      
    private:
 
      // Pseudofermion field phi
      LatticeFermion phi;

      // A handle for the EvenOddPrecWilsonFermAct
      Handle<const EvenOddPrecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;

      // The parameters for the inversion
      InvertParam_t inv_param;
    };


}; //end namespace chroma













#endif
