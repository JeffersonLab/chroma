// -*- C++ -*-
// $Id: central_tprec_logdet_tt_monomial_w.h,v 3.2 2008-05-23 18:39:45 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned log(det(A_ee))
 */

#ifndef __central_tprec_logdet_tt_monomial_h__
#define __central_tprec_logdet_tt_monomial_h__

#include "qdp_config.h"

#if QDP_NS == 4
#if QDP_NC == 3
#if QDP_ND == 4

#include "central_tprec_fermact_w.h"
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
  class CentralTimePrecLogDetTTMonomial : public ExactMonomial<P,Q>    
  {
  public: 
    virtual ~CentralTimePrecLogDetTTMonomial() {}

    void dsdq(P& F, const AbsFieldState<P,Q>& s) 
      {
	START_CODE();

	XMLWriter& xml_out = TheXMLLogWriter::Instance();
	push(xml_out, "CentralTimePrecLogDetTTMonomial");

	// Create FermAct
	const CentralTimePrecFermAct<Phi,P,Q>& FA = getFermAct();
      
	// Create a state for linop
	Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
	//Create LinOp
	Handle< CentralTimePrecLinearOperator<Phi,P,Q> > lin(FA.linOp(state));

	lin->derivLogDetTDagT(F, PLUS);
	F[lin->tDir()] *= Real(getNumFlavors());	  
	
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

	push(xml_out, "CentralTimePrecLogDetTTMonomial");
	const CentralTimePrecFermAct<Phi,P,Q>& FA = getFermAct();
	Handle< FermState<Phi,P,Q> > bc_g_state = FA.createState(s.getQ());

	// Need way to get gauge state from AbsFieldState<P,Q>
	Handle< CentralTimePrecLinearOperator<Phi,P,Q> > lin(FA.linOp(bc_g_state));
      
	Double S_ee =(Double(-2*getNumFlavors())*lin->logDetTDagT());
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
    virtual const CentralTimePrecFermAct<Phi,P,Q>& getFermAct() const = 0;
    virtual const int getNumFlavors() const = 0;
  };


  /*! @ingroup monomial */
  namespace CentralTimePrecLogDetTTMonomial4DEnv 
  {
    bool registerAll();
  }

  // Parameter structure
  /*! @ingroup monomial */
  struct CentralTimePrecLogDetTTMonomialParams 
  {
    // Base Constructor
    CentralTimePrecLogDetTTMonomialParams();

    // Read monomial from some root path
    CentralTimePrecLogDetTTMonomialParams(XMLReader& in, const std::string& path);
    GroupXML_t fermact;
    int num_flavors;
  };

  /*! @ingroup monomial */
  void read(XMLReader& r, const std::string& path,  CentralTimePrecLogDetTTMonomialParams& p);

  /*! @ingroup monomial */
  void write(XMLWriter& xml, const std::string& path, const CentralTimePrecLogDetTTMonomialParams& p);


  //! A Monomial For Just the EvenEven part of EvenOddPrecLogDetWilsonTypeFermActs -- concretely a 4D one
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class CentralTimePrecLogDetTTMonomial4D : 
    public CentralTimePrecLogDetTTMonomial<multi1d<LatticeColorMatrix>, 
    multi1d<LatticeColorMatrix>,
    LatticeFermion> 
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Construct from param struct
    CentralTimePrecLogDetTTMonomial4D(const CentralTimePrecLogDetTTMonomialParams& p);

    //! Copy Constructor
    CentralTimePrecLogDetTTMonomial4D(const CentralTimePrecLogDetTTMonomial4D& m) : num_flavors(m.num_flavors), fermact(m.fermact) {}

    //! Destructor is automagic
    ~CentralTimePrecLogDetTTMonomial4D() {}


  protected:
    const CentralTimePrecFermAct<T,P,Q> & getFermAct(void) const {
      return *fermact;
    }

    const int getNumFlavors() const {
      return num_flavors;
    }

  private:
    int num_flavors;
    Handle< CentralTimePrecFermAct<T,P,Q> > fermact;
  };

} //end namespace chroma

#endif
#endif
#endif

#endif
