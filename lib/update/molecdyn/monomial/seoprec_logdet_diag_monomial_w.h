// -*- C++ -*-
/*! \file
 *  \brief Symmetric even-odd preconditioned log(det(A_ee)) and log(det(A_oo))
 */

#ifndef __seoprec_logdet_diag_monomial_h__
#define __seoprec_logdet_diag_monomial_h__

#include "seoprec_logdet_wilstype_fermact_w.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/monomial/force_monitors.h"
#include "io/xmllog_io.h"


namespace Chroma 
{

  //! A Monomial For Just the diagonal parts of SymEvenOddPrecLogDetWilsonTypeFermActs
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  template<typename P, typename Q, typename Phi>
  class SymEvenOddPrecLogDetDiagMonomial : public ExactMonomial<P,Q>    
  {
  public: 
    virtual ~SymEvenOddPrecLogDetDiagMonomial() {}

    void dsdq(P& F, const AbsFieldState<P,Q>& s) 
    {
      START_CODE();

      XMLWriter& xml_out = TheXMLLogWriter::Instance();
      push(xml_out, "SymEvenOddPrecLogDetDiagMonomial");

      // Create FermAct
      const SymEvenOddPrecLogDetWilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();
      
      // Create a state for linop
      Handle< FermState<Phi,P,Q> > state(FA.createState(s.getQ()));
	
      //Create LinOp
      Handle< SymEvenOddPrecLogDetLinearOperator<Phi,P,Q> > lin(FA.linOp(state));

      // The derivative of each term
      P F_tmp;
      lin->derivLogDetEvenEvenLinOp(F, PLUS);
      lin->derivLogDetOddOddLinOp(F_tmp, PLUS);
      F += F_tmp;

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

      push(xml_out, "SymEvenOddPrecLogDetDiagMonomial");
      const SymEvenOddPrecLogDetWilsonTypeFermAct<Phi,P,Q>& FA = getFermAct();
      Handle< FermState<Phi,P,Q> > bc_g_state = FA.createState(s.getQ());

      // Need way to get gauge state from AbsFieldState<P,Q>
      Handle< SymEvenOddPrecLogDetLinearOperator<Phi,P,Q> > lin(FA.linOp(bc_g_state));
      
      Double S_ee = (Double(-getNumFlavors())*lin->logDetEvenEvenLinOp());
      Double S_oo = (Double(-getNumFlavors())*lin->logDetOddOddLinOp());
      Double S    = S_ee + S_oo;

      write(xml_out, "S", S);
      pop(xml_out);

      END_CODE();

      return S;
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
    virtual const SymEvenOddPrecLogDetWilsonTypeFermAct<Phi,P,Q>& getFermAct() const = 0;
    virtual const int getNumFlavors() const = 0;
  };


  /*! @ingroup monomial */
  namespace SymEvenOddPrecLogDetDiagMonomial4DEnv 
  {
    bool registerAll();
  }

  // Parameter structure
  /*! @ingroup monomial */
  struct SymEvenOddPrecLogDetDiagMonomialParams 
  {
    // Base Constructor
    SymEvenOddPrecLogDetDiagMonomialParams();

    // Read monomial from some root path
    SymEvenOddPrecLogDetDiagMonomialParams(XMLReader& in, const std::string& path);

    GroupXML_t fermact;
    int num_flavors;
  };

  /*! @ingroup monomial */
  void read(XMLReader& r, const std::string& path,  SymEvenOddPrecLogDetDiagMonomialParams& p);

  /*! @ingroup monomial */
  void write(XMLWriter& xml, const std::string& path, const SymEvenOddPrecLogDetDiagMonomialParams& p);


  //! A Monomial For Just the Diag part of SymEvenOddPrecLogDetWilsonTypeFermActs -- concretely a 4D one
  /*! @ingroup monomial
   *
   * Monomial is expected to be the same for these fermacts
   */
  class SymEvenOddPrecLogDetDiagMonomial4D : 
    public SymEvenOddPrecLogDetDiagMonomial<multi1d<LatticeColorMatrix>, 
					    multi1d<LatticeColorMatrix>,
					    LatticeFermion> 
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Construct from param struct
    SymEvenOddPrecLogDetDiagMonomial4D(const SymEvenOddPrecLogDetDiagMonomialParams& p);

    //! Copy Constructor
    SymEvenOddPrecLogDetDiagMonomial4D(const SymEvenOddPrecLogDetDiagMonomial4D& m) : num_flavors(m.num_flavors), fermact(m.fermact) {}

    //! Destructor is automagic
    ~SymEvenOddPrecLogDetDiagMonomial4D() {}


  protected:
    const SymEvenOddPrecLogDetWilsonTypeFermAct<T,P,Q> & getFermAct(void) const {
      return *fermact;
    }

    const int getNumFlavors() const {
      return num_flavors;
    }

  private:
    int num_flavors;
    Handle< SymEvenOddPrecLogDetWilsonTypeFermAct<T,P,Q> > fermact;
  };

} //end namespace chroma

#endif
