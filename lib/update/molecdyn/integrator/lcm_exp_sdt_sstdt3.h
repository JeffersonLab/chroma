// -*- C++ -*-

/*! @file
 * @brief Intgrator for exp(S dt - SST dt^3)
 *
 */

#ifndef LCM_EXPSDT_MINUS_SSTDT3_H
#define LCM_EXPSDT_MINUS_SSTDT3_H


#include "chromabase.h"
#include "handle.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/integrator/integrator_shared.h" 


namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatExpSdtMinusSSTdt3IntegratorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! @ingroup integrator */
  struct  LatColMatExpSdtMinusSSTdt3IntegratorParams
  {
    LatColMatExpSdtMinusSSTdt3IntegratorParams();
    LatColMatExpSdtMinusSSTdt3IntegratorParams(XMLReader& xml, const std::string& path);
    int  n_steps;
    multi1d<std::string> monomial_list;
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatExpSdtMinusSSTdt3IntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatExpSdtMinusSSTdt3IntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatExpSdtMinusSSTdt3Integrator 
    : public AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> > 
  {
  public:

    // Simplest Constructor
    LatColMatExpSdtMinusSSTdt3Integrator(const int  n_steps_, 
			      const multi1d<std::string>& monomial_id_list) : n_steps(n_steps_) { 
      create(monomial_id_list);
    }


    // Construct from nsteps and monomial_handle array
    LatColMatExpSdtMinusSSTdt3Integrator(const int n_steps_,
      const multi1d< Handle< Monomial< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > > monomials_) : n_steps(n_steps_), monomials(monomials_) {}

    // Construct from params struct and Hamiltonian
    LatColMatExpSdtMinusSSTdt3Integrator(const LatColMatExpSdtMinusSSTdt3IntegratorParams& p) : n_steps(p.n_steps) { 
      create(p.monomial_list);
    }

    // Copy constructor
    LatColMatExpSdtMinusSSTdt3Integrator(const LatColMatExpSdtMinusSSTdt3Integrator& l) :
      n_steps(l.n_steps), monomials(l.monomials) {}

    // ! Destruction is automagic
    ~LatColMatExpSdtMinusSSTdt3Integrator(void) {};

    void operator()( AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s, 
		     const Real& t1) const;

   			    

    void refreshFields(AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s) const { 
      for(int i=0; i < monomials.size(); i++) { 
	monomials[i]->refreshInternalFields(s);
      }
    }

    //! Reset Predictors in just this level
    void resetPredictors(void) const {
      for(int i=0; i < monomials.size(); ++i) {
	monomials[i]->resetPredictors();
      }
    }

  private:
    int  n_steps;

    typedef multi1d<LatticeColorMatrix> LCM;
    // Use the shared routine to bind the monomial list 
    void create(multi1d<std::string> monomial_id_list) { 
      IntegratorShared::bindMonomials(monomial_id_list, monomials);
    }

   

    multi1d< Handle< Monomial<multi1d<LatticeColorMatrix>, 
			      multi1d<LatticeColorMatrix> > > > monomials;
    
    
  };

}


#endif
