// -*- C++ -*-
// $Id: lcm_sts_force_grad_recursive.h,v 3.2 2017-12-20 11:20:00 boram Exp $

/*! @file
 * @brief Lat Col Mat force-gradient integrator
 *
 * Force-gradient integrator given in Kennedy and Clark, 2007 (arXiv:0710.3611)
 * Implementation adopted the Taylor expansion trick proposed in
 * Yin and Mawhinney, 2011 (arXiv:1111.5059)
 */

#ifndef LCM_FORCE_GRAD_RECURSIVE_H
#define LCM_FORCE_GRAD_RECURSIVE_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/integrator/integrator_shared.h"

namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatSTSForceGradRecursiveIntegratorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! @ingroup integrator */
  struct  LatColMatSTSForceGradRecursiveIntegratorParams
  {
    LatColMatSTSForceGradRecursiveIntegratorParams();
    LatColMatSTSForceGradRecursiveIntegratorParams(XMLReader& xml, const std::string& path);
    int  n_steps;
    multi1d<std::string> monomial_ids;
    std::string subintegrator_xml;
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatSTSForceGradRecursiveIntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatSTSForceGradRecursiveIntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatSTSForceGradRecursiveIntegrator 
    : public AbsRecursiveIntegrator<multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> > 
  {
  public:

    // Simplest Constructor
    LatColMatSTSForceGradRecursiveIntegrator(int  n_steps_, 
					 const multi1d<std::string>& monomial_ids_,
           
					 Handle< AbsComponentIntegrator< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& SubIntegrator_) : n_steps(n_steps_), SubIntegrator(SubIntegrator_) {

      IntegratorShared::bindMonomials(monomial_ids_, monomials);
    };

    // Construct from params struct and Hamiltonian
    LatColMatSTSForceGradRecursiveIntegrator(
					 const LatColMatSTSForceGradRecursiveIntegratorParams& p) : n_steps(p.n_steps), SubIntegrator(IntegratorShared::createSubIntegrator(p.subintegrator_xml)) {

      IntegratorShared::bindMonomials(p.monomial_ids, monomials);
      
    }


    // Copy constructor
    LatColMatSTSForceGradRecursiveIntegrator(const LatColMatSTSForceGradRecursiveIntegrator& l) :
      n_steps(l.n_steps), monomials(l.monomials), SubIntegrator(l.SubIntegrator) {}

    // ! Destruction is automagic
    ~LatColMatSTSForceGradRecursiveIntegrator(void) {};


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
    
    int  n_steps;

    multi1d< IntegratorShared::MonomialPair > monomials;

    Handle< AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
				   multi1d<LatticeColorMatrix> > > SubIntegrator;

    void fg_update( 
      AbsFieldState<multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >& s, 
      const Real& traj_length1, const Real& traj_length2) const;

  };

}


#endif
