// -*- C++ -*-
// $Id: lcm_exp_sdt.h,v 3.3 2006-12-28 17:34:00 bjoo Exp $

/*! @file
 * @brief Intgrator for exp(S dt)
 *
 * A component integrator to integrate exp(Sdt) for some 
 * monomial list S (ie this is a leapP like component)
 */

#ifndef LCM_EXPSDT_H
#define LCM_EXPSDT_H


#include "chromabase.h"
#include "handle.h"
#include "update/molecdyn/monomial/abs_monomial.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/integrator/integrator_shared.h" 


namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatExpSdtIntegratorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! @ingroup integrator */
  struct  LatColMatExpSdtIntegratorParams
  {
    LatColMatExpSdtIntegratorParams();
    LatColMatExpSdtIntegratorParams(XMLReader& xml, const std::string& path);
    int  n_steps;
    multi1d<std::string> monomial_list;
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatExpSdtIntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatExpSdtIntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatExpSdtIntegrator 
    : public AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> > 
  {
  public:

    // Simplest Constructor
    LatColMatExpSdtIntegrator(const int  n_steps_, 
			      const multi1d<std::string>& monomial_id_list) : n_steps(n_steps_) { 
      create(monomial_id_list);
    }


    // Construct from nsteps and monomial_handle array
    LatColMatExpSdtIntegrator(const int n_steps_,
			      const multi1d< IntegratorShared::MonomialPair>& monomials_) : n_steps(n_steps_), monomials(monomials_) {}

    // Construct from params struct and Hamiltonian
    LatColMatExpSdtIntegrator(const LatColMatExpSdtIntegratorParams& p) : n_steps(p.n_steps) { 
      create(p.monomial_list);
    }

    // Copy constructor
    LatColMatExpSdtIntegrator(const LatColMatExpSdtIntegrator& l) :
      n_steps(l.n_steps), monomials(l.monomials) {}

    // ! Destruction is automagic
    ~LatColMatExpSdtIntegrator(void) {};

    void operator()( AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s, 
		     const Real& traj_length) const;
   			    

    void refreshFields(AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s) const { 
      for(int i=0; i < monomials.size(); i++) { 
	monomials[i].mon->refreshInternalFields(s);
      }
    }

    //! Reset Predictors in just this level
    void resetPredictors(void) const {
      for(int i=0; i < monomials.size(); ++i) {
	monomials[i].mon->resetPredictors();
      }
    }

  private:
    int  n_steps;

    typedef multi1d<LatticeColorMatrix> LCM;
    // Use the shared routine to bind the monomial list 
    void create(const multi1d<std::string>& monomial_id_list) { 
      IntegratorShared::bindMonomials(monomial_id_list, monomials);
    }

   

    multi1d< IntegratorShared::MonomialPair > monomials;
  };

}


#endif
