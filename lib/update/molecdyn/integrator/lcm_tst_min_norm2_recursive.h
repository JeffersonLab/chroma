// -*- C++ -*-
// $Id: lcm_tst_min_norm2_recursive.h,v 3.2 2006-12-28 17:34:00 bjoo Exp $

/*! @file
 * @brief Intgrator for exp(S dt)
 *
 * A component integrator to integrate exp(Sdt) for some 
 * monomial list S (ie this is a leapP like component)
 */

#ifndef LCM_TST_MIN_NORM2_RECURSIVE_H
#define LCM_TST_MIN_NORM2_RECURSIVE_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/integrator/integrator_shared.h"

namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatTSTMinNorm2RecursiveIntegratorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! @ingroup integrator */
  struct  LatColMatTSTMinNorm2RecursiveIntegratorParams
  {
    LatColMatTSTMinNorm2RecursiveIntegratorParams();
    LatColMatTSTMinNorm2RecursiveIntegratorParams(XMLReader& xml, const std::string& path);
    int  n_steps;
    Real lambda;
    multi1d<std::string> monomial_ids;
    std::string subintegrator_xml;
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatTSTMinNorm2RecursiveIntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatTSTMinNorm2RecursiveIntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatTSTMinNorm2RecursiveIntegrator 
    : public AbsRecursiveIntegrator<multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> > 
  {
  public:

    // Simplest Constructor
    LatColMatTSTMinNorm2RecursiveIntegrator(int  n_steps_, 
					 const multi1d<std::string>& monomial_ids_,
					 Real lambda_, 

					 Handle< AbsComponentIntegrator< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& SubIntegrator_) : n_steps(n_steps_), lambda(lambda_), SubIntegrator(SubIntegrator_) {

      IntegratorShared::bindMonomials(monomial_ids_, monomials);
    };

    // Construct from params struct and Hamiltonian
    LatColMatTSTMinNorm2RecursiveIntegrator(
					 const LatColMatTSTMinNorm2RecursiveIntegratorParams& p) : n_steps(p.n_steps), lambda(p.lambda), SubIntegrator(IntegratorShared::createSubIntegrator(p.subintegrator_xml)) {

      IntegratorShared::bindMonomials(p.monomial_ids, monomials);
      
    }


    // Copy constructor
    LatColMatTSTMinNorm2RecursiveIntegrator(const LatColMatTSTMinNorm2RecursiveIntegrator& l) :
      n_steps(l.n_steps), monomials(l.monomials), lambda(l.lambda), SubIntegrator(l.SubIntegrator) {}

    // ! Destruction is automagic
    ~LatColMatTSTMinNorm2RecursiveIntegrator(void) {};


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
    Real lambda;

    multi1d< IntegratorShared::MonomialPair > monomials;

    Handle< AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
				   multi1d<LatticeColorMatrix> > > SubIntegrator;
	              

  };

}


#endif
