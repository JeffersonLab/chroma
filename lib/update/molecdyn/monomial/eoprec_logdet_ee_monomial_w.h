// -*- C++ -*-
// $Id: eoprec_logdet_ee_monomial_w.h,v 3.4 2008-05-23 18:39:45 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned log(det(A_ee))
 */

#ifndef __eoprec_logdet_even_even_monomial_h__
#define __eoprec_logdet_even_even_monomial_h__

#include "eoprec_logdet_wilstype_fermact_w.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "io/xmllog_io.h"


namespace Chroma 
{

  //! A Monomial For Just the EvenEven part of EvenOddPrecLogDetWilsonTypeFermActs
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  template<typename P, typename Q, typename Phi>
  class EvenOddPrecLogDetEvenEvenMonomial : public ExactMonomial<P,Q>    
  {
  public: 
    virtual ~EvenOddPrecLogDetEvenEvenMonomial() {}

    void dsdq(P& F, const AbsFieldState<P,Q>& s) 
      {
	START_CODE();

	XMLWriter& xml_out = TheXMLLogWriter::Instance();
	push(xml_out, "EvenOddPrecLogDetEvenEvenMonomial");

	// Create FermAct
	const EvenOddPrecLogDetWilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();
      
	// Create a state for linop
	Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
	//Create LinOp
	Handle< EvenOddPrecLogDetLinearOperator<Phi,P,Q> > lin(FA.linOp(state));

	lin->derivLogDetEvenEvenLinOp(F, PLUS);

	for(int mu=0; mu < Nd; mu++) { 
	  F[mu] *= Real(-getNumFlavors());	  
	}
	
	state->deriv(F);

	monitorForces(xml_out, "Forces", F);
	pop(xml_out);

	END_CODE();
      }


    //! Gauge action value
    Double S(const AbsFieldState<P,Q>& s)  
      {
	START_CODE();

	XMLWriter& xml_out = TheXMLLogWriter::Instance();

	push(xml_out, "EvenOddPrecLogDetEvenEvenMonomial");
	const EvenOddPrecLogDetWilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();
	Handle< FermState<Phi,P,Q> > bc_g_state = FA.createState(s.getQ());

	// Need way to get gauge state from AbsFieldState<P,Q>
	Handle< EvenOddPrecLogDetLinearOperator<Phi,P,Q> > lin(FA.linOp(bc_g_state));
      
	Double S_ee =(Double(-getNumFlavors())*lin->logDetEvenEvenLinOp());
	write(xml_out, "S", S_ee);
	pop(xml_out);

	END_CODE();

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
  namespace EvenOddPrecLogDetEvenEvenMonomial4DEnv 
  {
    bool registerAll();
  }

  // Parameter structure
  /*! @ingroup monomial */
  struct EvenOddPrecLogDetEvenEvenMonomialParams 
  {
    // Base Constructor
    EvenOddPrecLogDetEvenEvenMonomialParams();

    // Read monomial from some root path
    EvenOddPrecLogDetEvenEvenMonomialParams(XMLReader& in, const std::string& path);
    GroupXML_t fermact;
    int num_flavors;
  };

  /*! @ingroup monomial */
  void read(XMLReader& r, const std::string& path,  EvenOddPrecLogDetEvenEvenMonomialParams& p);

  /*! @ingroup monomial */
  void write(XMLWriter& xml, const std::string& path, const EvenOddPrecLogDetEvenEvenMonomialParams& p);


  //! A Monomial For Just the EvenEven part of EvenOddPrecLogDetWilsonTypeFermActs -- concretely a 4D one
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class EvenOddPrecLogDetEvenEvenMonomial4D : 
    public EvenOddPrecLogDetEvenEvenMonomial<multi1d<LatticeColorMatrix>, 
    multi1d<LatticeColorMatrix>,
    LatticeFermion> 
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Construct from param struct
    EvenOddPrecLogDetEvenEvenMonomial4D(const EvenOddPrecLogDetEvenEvenMonomialParams& p);

    //! Copy Constructor
    EvenOddPrecLogDetEvenEvenMonomial4D(const EvenOddPrecLogDetEvenEvenMonomial4D& m) : num_flavors(m.num_flavors), fermact(m.fermact) {}

    //! Destructor is automagic
    ~EvenOddPrecLogDetEvenEvenMonomial4D() {}


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

} //end namespace chroma

#endif
