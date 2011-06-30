// -*- C++ -*-
// $Id: lcm_tst_min_norm2_recursive_dtau.h,v 3.1 2008-11-16 17:04:04 bjoo Exp $

/*! @file
 * @brief Intgrator for exp(S dt)
 *
 * A component integrator to integrate exp(Sdt) for some 
 * monomial list S (ie this is a leapP like component)
 */

#ifndef LCM_TST_MIN_NORM2_RECURSIVE_DTAU_H
#define LCM_TST_MIN_NORM2_RECURSIVE_DTAU_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/integrator/integrator_shared.h"

namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatTSTMinNorm2DTauRecursiveIntegratorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! @ingroup integrator */
  struct  LatColMatTSTMinNorm2DTauRecursiveIntegratorParams
  {
    LatColMatTSTMinNorm2DTauRecursiveIntegratorParams();
    LatColMatTSTMinNorm2DTauRecursiveIntegratorParams(XMLReader& xml, const std::string& path);
    Real delta_tau_max; // Desired delta_tau_max
    Real lambda;
    multi1d<std::string> monomial_ids;
    std::string subintegrator_xml;
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatTSTMinNorm2DTauRecursiveIntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatTSTMinNorm2DTauRecursiveIntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatTSTMinNorm2DTauRecursiveIntegrator 
    : public AbsRecursiveIntegrator<multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> > 
  {
  public:

    // Simplest Constructor
    LatColMatTSTMinNorm2DTauRecursiveIntegrator(Real  delta_tau_max_, 
					 const multi1d<std::string>& monomial_ids_,
					 Real lambda_, 

					 Handle< AbsComponentIntegrator< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& SubIntegrator_) : delta_tau_max(delta_tau_max_), lambda(lambda_), SubIntegrator(SubIntegrator_) {

      IntegratorShared::bindMonomials(monomial_ids_, monomials);
    };

    // Construct from params struct and Hamiltonian
    LatColMatTSTMinNorm2DTauRecursiveIntegrator(
					 const LatColMatTSTMinNorm2DTauRecursiveIntegratorParams& p) : delta_tau_max(p.delta_tau_max), lambda(p.lambda), SubIntegrator(IntegratorShared::createSubIntegrator(p.subintegrator_xml)) {

      IntegratorShared::bindMonomials(p.monomial_ids, monomials);
      
    }


    // Copy constructor
    LatColMatTSTMinNorm2DTauRecursiveIntegrator(const LatColMatTSTMinNorm2DTauRecursiveIntegrator& l) :
      delta_tau_max(l.delta_tau_max), monomials(l.monomials), lambda(l.lambda), SubIntegrator(l.SubIntegrator) {}

    // ! Destruction is automagic
    ~LatColMatTSTMinNorm2DTauRecursiveIntegrator(void) {};


    void operator()( AbsFieldState<multi1d<LatticeColorMatrix>,
		                   multi1d<LatticeColorMatrix> >& s, 
		     const Real& traj_length) const;
   			    
    AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
			   multi1d<LatticeColorMatrix> >& getSubIntegrator() const {
      return (*SubIntegrator);
    }
    
  protected:
    //! Refresh fields in just this level
    void refreshFieldsThisLevel(AbsFieldState<multi1d<LatticeColorMatrix>,
				multi1d<LatticeColorMatrix> >& s) const {
      for(int i=0; i < monomials.size(); i++) { 
	monomials[i].mon->refreshInternalFields(s);
      }
    }

    //! Reset Predictors in just this level
    void resetPredictorsThisLevel(void) const {
      for(int i=0; i < monomials.size(); ++i) {
	monomials[i].mon->resetPredictors();
      }
    }

  private:
    
  Real delta_tau_max;  // Maximum delta_tau (may be shortened to make n_steps fit
                   // but this is the maximum desired timescale
  Real lambda;
  
    multi1d< IntegratorShared::MonomialPair > monomials;
    
    Handle< AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
				   multi1d<LatticeColorMatrix> > > SubIntegrator;
	              

  };

}


#endif
