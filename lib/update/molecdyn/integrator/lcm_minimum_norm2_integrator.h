
/*! @file
 * @brief 2nd order minimum norm intergrator a la 
 * Omelyan adapted to QCD by deForcrand and Takaishi
 * hep-lat/0505020
 *
 * Leapfrog intregator for HMC
 */

#ifndef LCM_MINIMUM_NORM_H
#define LCM_MINIMUM_NORM_H


#include "chromabase.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/integrator/lcm_integrator_leaps.h"
#include <string>


namespace Chroma 
{

  /*! @ingroup integrator */
  namespace LatColMatMinimumNorm2IntegratorEnv {
    extern const std::string name;
    extern const bool registered;
  };

  /*! @ingroup integrator */
  struct  LatColMatMinimumNorm2IntegratorParams {
    LatColMatMinimumNorm2IntegratorParams();
    LatColMatMinimumNorm2IntegratorParams(XMLReader& xml, const std::string& path);

    int n_steps;
    Real tau0;
    Real lambda; 

  };

  /*! @ingroup integrator */
  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatMinimumNorm2IntegratorParams& p);

  /*! @ingroup integrator */
  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatMinimumNorm2IntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  /*! @ingroup integrator
   *  Specialised to multi1d<LatticeColorMatrix>
   */
  class LatColMatMinimumNorm2Integrator 
    : public AbsMDIntegrator<multi1d<LatticeColorMatrix>,
                             multi1d<LatticeColorMatrix> > {

  public:

    // Construct from params struct and Hamiltonian
    LatColMatMinimumNorm2Integrator(const LatColMatMinimumNorm2IntegratorParams& p,
				    Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_) : n_steps(p.n_steps), lambda(p.lambda), tau0(p.tau0), H_MD(H_) { 

    }

    // Copy constructor
    LatColMatMinimumNorm2Integrator(const LatColMatMinimumNorm2Integrator& l) :
      n_steps(l.n_steps), lambda(l.lambda), tau0(l.tau0), H_MD(l.H_MD) {}

    // ! Destruction is automagic
    ~LatColMatMinimumNorm2Integrator(void) {};

    //! Get at the MD Hamiltonian
    AbsHamiltonian<multi1d<LatticeColorMatrix>,
		   multi1d<LatticeColorMatrix> >& getHamiltonian(void) {
      return *H_MD;
    }

    //! Do a trajectory
    void operator()(AbsFieldState<multi1d<LatticeColorMatrix>,
		    multi1d<LatticeColorMatrix> >& s) {

      Real dt = getStepSize();      // Overall step size
      Real dtby2 = dt/Real(2);
      Real lambdadt = lambda*dt;
      Real twolambdadt = Real(2)*lambdadt;
      Real one_m_twolambda_dt=(Real(1)-Real(2)*lambda)*dt;
      
      leapP(lambdadt,s);
      leapQ(dtby2,s);
      leapP(one_m_twolambda_dt, s);
      leapQ(dtby2,s);
      for(int step=0; step < n_steps-1; step++) {
		leapP(twolambdadt,s);
		leapQ(dtby2,s);
		leapP(one_m_twolambda_dt, s);
		leapQ(dtby2,s);
      }
      leapP(lambdadt,s);

    }

  protected:

    //! LeapP for just a selected list of monomials
    void leapP(const Real& dt, 
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	       multi1d<LatticeColorMatrix> >& s) {
      LCMMDIntegratorSteps::leapP(dt, getHamiltonian(), s);
    }


    //! Leap with Q (with all monomials)
    void leapQ(const Real& dt, 
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	       multi1d<LatticeColorMatrix> >& s) {
      LCMMDIntegratorSteps::leapQ(dt, s);
    }


    //! Get the trajectory length
    const Real getTrajLength(void) const {
      return tau0;
    }

    //! Get the step size 
    const Real getStepSize(void) const {
      return tau0/Real(n_steps);
    }

    
  private:
    int n_steps;
    Real lambda;
    Real tau0;
    Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, 
			    multi1d<LatticeColorMatrix> > > H_MD;
    
  };

};


#endif
