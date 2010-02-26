// -*- C++ -*-


#ifndef LCM_FORCE_GRADIENT_RECURSIVE_H
#define LCM_FORCE_GRADIENT_RECURSIVE_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/integrator/integrator_shared.h"

namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatSTSForceGradientRecursiveIntegratorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  /*! @ingroup integrator */
  struct  LatColMatSTSForceGradientRecursiveIntegratorParams
  {
    LatColMatSTSForceGradientRecursiveIntegratorParams();
    LatColMatSTSForceGradientRecursiveIntegratorParams(XMLReader& xml, const std::string& path);
    int  n_steps;
    multi1d<std::string> monomial_ids;
    std::string subintegrator_xml;
  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatSTSForceGradientRecursiveIntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatSTSForceGradientRecursiveIntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatSTSForceGradientRecursiveIntegrator 
    : public AbsRecursiveIntegrator<multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> > 
  {
  public:

    // Simplest Constructor
    LatColMatSTSForceGradientRecursiveIntegrator(int  n_steps_, 
					 const multi1d<std::string>& monomial_ids_,

					 Handle< AbsComponentIntegrator< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& SubIntegrator_) : n_steps(n_steps_), SubIntegrator(SubIntegrator_) {

      IntegratorShared::bindMonomials(monomial_ids_, monomials);
    };

    // Construct from params struct and Hamiltonian
    LatColMatSTSForceGradientRecursiveIntegrator(
					 const LatColMatSTSForceGradientRecursiveIntegratorParams& p) : n_steps(p.n_steps), SubIntegrator(IntegratorShared::createSubIntegrator(p.subintegrator_xml)) {

      IntegratorShared::bindMonomials(p.monomial_ids, monomials);
      
    }


    // Copy constructor
    LatColMatSTSForceGradientRecursiveIntegrator(const LatColMatSTSForceGradientRecursiveIntegrator& l) :
      n_steps(l.n_steps), monomials(l.monomials), SubIntegrator(l.SubIntegrator) {}

    // ! Destruction is automagic
    ~LatColMatSTSForceGradientRecursiveIntegrator(void) {};


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

    //! Reset Predictors in just this level
    void resetPredictorsThisLevel(void) const {
      for(int i=0; i < monomials.size(); ++i) {
	monomials[i]->resetPredictors();
      }
    }

  private:
    
    int  n_steps;

    multi1d< Handle< Monomial< multi1d<LatticeColorMatrix>, 
			       multi1d<LatticeColorMatrix> > > > monomials;

    Handle< AbsComponentIntegrator<multi1d<LatticeColorMatrix>,
				   multi1d<LatticeColorMatrix> > > SubIntegrator;
	              
    void getShadow( AbsFieldState<multi1d<LatticeColorMatrix>,
		    multi1d<LatticeColorMatrix> >& s,
		    const Real& dt,
		    Double& H, 
		    Double& Hs) const;
  };

}


#endif
