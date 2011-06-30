// -*- C++ -*-
// $Id: lcm_creutz_gocksch_4_recursive.h,v 1.2 2006-12-28 17:34:00 bjoo Exp $

/*! @file
 * @brief Lat Col Mat 4th order Creutz-Gocksch (Campostrini?) Integrator
 *
 * Integrator that is accurate to 5th order per time step or 4th order 
 * per trajectory. Creutz-Gocksch (Campostrini?) Construction. This
 * particular one is from the paper of Sexton and Weingarten eq: 6.1
 */

#ifndef LCM_CREUTZ_GOCKSCH_RECURSIVE_4_H
#define LCM_CREUTZ_GOCKSCH_RECURSIVE_4_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/integrator/integrator_shared.h"

namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatCreutzGocksch4RecursiveIntegratorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! @ingroup integrator */
  struct  LatColMatCreutzGocksch4RecursiveIntegratorParams
  {
    LatColMatCreutzGocksch4RecursiveIntegratorParams();
    LatColMatCreutzGocksch4RecursiveIntegratorParams(XMLReader& xml, const std::string& path);
    int  n_steps;

    multi1d<std::string> monomial_ids;
    std::string subintegrator_xml;
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatCreutzGocksch4RecursiveIntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatCreutzGocksch4RecursiveIntegratorParams& p);

  //! MD integrator interface for 4th order 5 Force Min. Norm. Integrator (Velocity variant)
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatCreutzGocksch4RecursiveIntegrator 
    : public AbsRecursiveIntegrator<multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> > 
  {
  public:
    
    // Construct from params struct and Hamiltonian
    LatColMatCreutzGocksch4RecursiveIntegrator(
               const LatColMatCreutzGocksch4RecursiveIntegratorParams& p
	  )  : params(p), SubIntegrator(IntegratorShared::createSubIntegrator(p.subintegrator_xml)) {
      IntegratorShared::bindMonomials(p.monomial_ids, monomials);
    }


    // Copy constructor
    LatColMatCreutzGocksch4RecursiveIntegrator(const LatColMatCreutzGocksch4RecursiveIntegrator& l) :
      params(l.params), monomials(l.monomials), SubIntegrator(l.SubIntegrator) {}

    // ! Destruction is automagic
    ~LatColMatCreutzGocksch4RecursiveIntegrator(void) {};


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
    LatColMatCreutzGocksch4RecursiveIntegratorParams params;


  multi1d< IntegratorShared::MonomialPair > monomials;

    Handle< AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
				   multi1d<LatticeColorMatrix> > > SubIntegrator;
	              

    
  };

}


#endif
