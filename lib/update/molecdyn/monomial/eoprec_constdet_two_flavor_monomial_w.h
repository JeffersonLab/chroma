// -*- C++ -*-
// $Id: eoprec_constdet_two_flavor_monomial_w.h,v 3.2 2008-05-21 17:07:50 bjoo Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#ifndef __prec_two_flavor_monomial_w_h__
#define __prec_two_flavor_monomial_w_h__

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/two_flavor_monomial_w.h"
#include "update/molecdyn/monomial/two_flavor_monomial_params_w.h"

namespace Chroma 
{

  /*! @ingroup monomial */
  namespace EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomialEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Wrapper class for  2-flavor even-odd prec ferm monomials
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial :
    public  TwoFlavorExactEvenOddPrecConstDetWilsonTypeFermMonomial< 
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix>,
    LatticeFermion>
  {
  public: 
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    // Construct out of a parameter struct. Check against the desired FermAct name
    EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial(const TwoFlavorWilsonTypeFermMonomialParams& param_);

    // Copy Constructor
    EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial(const EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial& m) : phi(m.phi), fermact(m.fermact), inv_param(m.inv_param), chrono_predictor(m.chrono_predictor) {}


    //! Even even contribution (eg ln det Clover)
    Double S_even_even(const AbsFieldState<multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >& s) {
      return Double(0);
    }

#if 1
    //! Compute the odd odd contribution (eg
    Double S_odd_odd(const AbsFieldState<multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >& s) {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "S_odd_odd");

      const EvenOddPrecWilsonTypeFermAct<T,P,Q>& FA = getFermAct();

      Handle< FermState<T,P,Q> > bc_g_state = FA.createState(s.getQ());

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< EvenOddPrecLinearOperator<T,P,Q> > lin(FA.linOp(bc_g_state));
      // Get the X fields
      T X;

      // Action calc doesnt use chrono predictor use zero guess
      X[ lin->subset() ] = zero;

      // getX noe always uses chrono predictor. Best to Nuke it therefore
      QDPIO::cout << "TwoFlavWilson4DMonomial: resetting Predictor before energy calc solve" << endl;
      (getMDSolutionPredictor()).reset();
      int n_count = getX(X, s);

      LatticeDouble site_action=zero;
      site_action[ lin->subset() ] = Double(-12);
      site_action[ lin->subset() ] += localInnerProductReal(getPhi(),X);

     
      //Double action = innerProductReal(getPhi(), X, lin->subset());
      Double action = sum(site_action, lin->subset());

      write(xml_out, "n_count", n_count);
      write(xml_out, "S_oo", action);
      pop(xml_out);

      END_CODE();

      return action;
    }
#endif
  protected:

    T& getPhi(void) {
      return phi;
    }

    const T& getPhi(void) const {
      return phi;
    }

    const EvenOddPrecWilsonTypeFermAct<T,P,Q>& getFermAct(void) const { 
      return *fermact;
    }
      
    //! Get parameters for the inverter
    const GroupXML_t& getInvParams(void) const { 
      return inv_param;
    }

    AbsChronologicalPredictor4D<T>& getMDSolutionPredictor(void) { 
      return *chrono_predictor;
    };


  private:
 
    // Hide empty constructor and =
    EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial();
    void operator=(const EvenOddPrecConstDetTwoFlavorWilsonTypeFermMonomial&);

    // Pseudofermion field phi
    T phi;

    // A handle for the EvenOddPrecWilsonFermAct
    Handle<const EvenOddPrecWilsonTypeFermAct<T,P,Q> > fermact;

    // The parameters for the inversion
    GroupXML_t inv_param;

    Handle<AbsChronologicalPredictor4D<T> > chrono_predictor;
    };


} //end namespace chroma

#endif
