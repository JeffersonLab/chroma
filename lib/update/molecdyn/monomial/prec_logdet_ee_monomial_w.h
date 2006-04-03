// -*- C++ -*-
// $Id: prec_logdet_ee_monomial_w.h,v 3.0 2006-04-03 04:59:09 edwards Exp $
/*! \file
 *  \brief Generic gauge action monomial wrapper
 */

#ifndef __prec_logdet_even_even_monomial_h__
#define __prec_logdet_even_even_monomial_h__

#include "chromabase.h"

#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/abs_monomial.h"

#include "io/xmllog_io.h"

namespace Chroma 
{

  //! A Monomial For Just the EvenEven part of EvenOddPrecLogDetWilsonTypeFermActs
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  template<typename P, typename Q, typename Phi>
  class PrecLogDetEvenEvenMonomial : public ExactMonomial<P,Q>    
  {
  public: 
    virtual ~PrecLogDetEvenEvenMonomial() {}

    void dsdq(P& F, const AbsFieldState<P,Q>& s) 
      {
	XMLWriter& xml_out = TheXMLOutputWriter::Instance();
	push(xml_out, "PrecLogDetEvenEvenMonomial");

	// Create FermAct
	const EvenOddPrecLogDetWilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();
      
	// Create a state for linop
	Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
	//Create LinOp
	Handle< EvenOddPrecLogDetLinearOperator<Phi,P,Q> > lin(FA.linOp(state));

	lin->derivEvenEvenLogDet(F, PLUS);

	for(int mu=0; mu < Nd; mu++) { 
	  F[mu] *= Real(-getNumFlavors());	  
	}
	
	state->deriv(F);

	Double F_sq = norm2(F);
	write(xml_out, "F_sq", F_sq);
	pop(xml_out);
      }


    //! Gauge action value
    Double S(const AbsFieldState<P,Q>& s)  {

      XMLWriter& xml_out = TheXMLOutputWriter::Instance();

      push(xml_out, "PrecLogDetEvenEvenMonomial");
      const EvenOddPrecLogDetWilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();
      Handle< FermState<Phi,P,Q> > bc_g_state = FA.createState(s.getQ());

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< EvenOddPrecLogDetLinearOperator<Phi,P,Q> > lin(FA.linOp(bc_g_state));
      
      Double S_ee =(Double(-getNumFlavors())*lin->LogDetEvenEven());
      write(xml_out, "S", S_ee);
      pop(xml_out);

      return S_ee;
    }
	
	
    void refreshInternalFields(const AbsFieldState<multi1d<LatticeColorMatrix>,
			       multi1d<LatticeColorMatrix> >& s) {
      //No internal fields to refresh => Nop
    }

    void setInternalFields(const Monomial<multi1d<LatticeColorMatrix>, 
			   multi1d<LatticeColorMatrix> >& m) {
      // No internal fields to refresh => Nop
    }

  protected:
    virtual const EvenOddPrecLogDetWilsonTypeFermAct<Phi,P,Q>& getFermAct() const = 0;
    virtual const int getNumFlavors() const = 0;
  };


  /*! @ingroup monomial */
  namespace PrecLogDetEvenEvenMonomial4DEnv {
    extern const std::string name;
    extern const bool registered;
  };

  // Parameter structure
  /*! @ingroup monomial */
  struct PrecLogDetEvenEvenMonomialParams 
  {
    // Base Constructor
    PrecLogDetEvenEvenMonomialParams();

    // Read monomial from some root path
    PrecLogDetEvenEvenMonomialParams(XMLReader& in, const std::string& path);
    std::string ferm_act;
    int num_flavors;
  };

  void read(XMLReader& r, const std::string& path,  PrecLogDetEvenEvenMonomialParams& p);

  void write(XMLWriter& xml, const std::string& path, const PrecLogDetEvenEvenMonomialParams& p);


  //! A Monomial For Just the EvenEven part of EvenOddPrecLogDetWilsonTypeFermActs -- concretely a 4D one
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class PrecLogDetEvenEvenMonomial4D : 
    public PrecLogDetEvenEvenMonomial<multi1d<LatticeColorMatrix>, 
                                      multi1d<LatticeColorMatrix>,
                                      LatticeFermion> 
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Construct from param struct
    PrecLogDetEvenEvenMonomial4D(const PrecLogDetEvenEvenMonomialParams& p);

    //! Copy Constructor
    PrecLogDetEvenEvenMonomial4D(const PrecLogDetEvenEvenMonomial4D& m) : num_flavors(m.num_flavors), fermact(m.fermact) {}

    //! Destructor is automagic
    ~PrecLogDetEvenEvenMonomial4D() {}


  protected:
    const EvenOddPrecLogDetWilsonTypeFermAct<T,P,Q> & getFermAct(void) const {
      return *fermact;
    }

    const int getNumFlavors() const {
      return num_flavors;
    }

  private:
    int num_flavors;
    Handle< EvenOddPrecLogDetWilsonTypeFermAct<T,P,Q> > fermact;
  };
}; //end namespace chroma

#endif
