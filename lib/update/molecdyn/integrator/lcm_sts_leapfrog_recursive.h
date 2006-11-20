// -*- C++ -*-
// $Id: lcm_sts_leapfrog_recursive.h,v 3.2 2006-11-20 19:15:02 bjoo Exp $

/*! @file
 * @brief Intgrator for exp(S dt)
 *
 * A component integrator to integrate exp(Sdt) for some 
 * monomial list S (ie this is a leapP like component)
 */

#ifndef LCM_STS_LEAPFROG_RECURSIVE_H
#define LCM_STS_LEAPFROG_RECURSIVE_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/integrator/integrator_shared.h"

namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatSTSLeapfrogRecursiveIntegratorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! @ingroup integrator */
  struct  LatColMatSTSLeapfrogRecursiveIntegratorParams
  {
    LatColMatSTSLeapfrogRecursiveIntegratorParams();
    LatColMatSTSLeapfrogRecursiveIntegratorParams(XMLReader& xml, const std::string& path);
    int  n_steps;
    multi1d<std::string> monomial_ids;
    std::string subintegrator_xml;
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatSTSLeapfrogRecursiveIntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatSTSLeapfrogRecursiveIntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatSTSLeapfrogRecursiveIntegrator 
    : public AbsRecursiveIntegrator<multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> > 
  {
  public:
    
    // Simplest Constructor
    LatColMatSTSLeapfrogRecursiveIntegrator(int  n_steps_, 
					    const multi1d<std::string>& monomial_ids_,
					    Handle< AbsComponentIntegrator< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& SubIntegrator_) : n_steps(n_steps_), SubIntegrator(SubIntegrator_) {
      IntegratorShared::bindMonomials(monomial_ids_, monomials);
    };

    // Construct from params struct and Hamiltonian
    LatColMatSTSLeapfrogRecursiveIntegrator(
               const LatColMatSTSLeapfrogRecursiveIntegratorParams& p
	  )  : n_steps(p.n_steps), SubIntegrator(IntegratorShared::createSubIntegrator(p.subintegrator_xml)) {
      IntegratorShared::bindMonomials(p.monomial_ids, monomials);
    }


    // Copy constructor
    LatColMatSTSLeapfrogRecursiveIntegrator(const LatColMatSTSLeapfrogRecursiveIntegrator& l) :
      n_steps(l.n_steps), monomials(l.monomials), SubIntegrator(l.SubIntegrator) {}

    // ! Destruction is automagic
    ~LatColMatSTSLeapfrogRecursiveIntegrator(void) {};


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
	monomials[i]->refreshInternalFields(s);
      }
    }
  private:

    int  n_steps;

    multi1d< Handle< Monomial< multi1d<LatticeColorMatrix>,
			       multi1d<LatticeColorMatrix> > > > monomials;

    Handle< AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
				   multi1d<LatticeColorMatrix> > > SubIntegrator;
	              

    
  };

}


#endif
