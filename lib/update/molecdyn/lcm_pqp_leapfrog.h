#ifndef LCM_PQP_LEAPFROG_H
#define LCM_PQP_LEAPFROG_H


#include "chromabase.h"
#include "abs_hamiltonian.h"
#include "abs_integrator.h"

#include <string>

using namespace QDP;
using namespace Chroma;
using namespace std;

namespace Chroma {

  namespace LatColMatPQPLeapfrogIntegratorEnv {
    extern const std::string name;
    extern const bool registered;
  };

  struct  LatColMatPQPLeapfrogIntegratorParams {
    LatColMatPQPLeapfrogIntegratorParams();
    LatColMatPQPLeapfrogIntegratorParams(XMLReader& xml, const std::string& path);

    Real dt;
    Real tau0;
  };

  void read(XMLReader& xml_in, 
	    const std::string& path,
	    LatColMatPQPLeapfrogIntegratorParams& p);

  void write(XMLWriter& xml_out,
	     const std::string& path, 
	     const LatColMatPQPLeapfrogIntegratorParams& p);

  //! MD integrator interface for PQP leapfrog
  //  specialised to multi1d<LatticeColorMatrix>
  class LatColMatPQPLeapfrogIntegrator 
    : public PQPLeapfrogIntegrator<multi1d<LatticeColorMatrix>,
				   multi1d<LatticeColorMatrix> > {

  public:

    // Simplest Constructor
    LatColMatPQPLeapfrogIntegrator(const Real& dt_, 
				   const Real& tau_,
				   Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_ ) : dt(dt_), tau0(tau_), H_MD(H_) {};

    // Construct from params struct and Hamiltonian
    LatColMatPQPLeapfrogIntegrator(const LatColMatPQPLeapfrogIntegratorParams& p,
				   Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& H_) : dt(p.dt), tau0(p.tau0), H_MD(H_) {}

    // Copy constructor
    LatColMatPQPLeapfrogIntegrator(const LatColMatPQPLeapfrogIntegrator& l) :
      dt(l.dt), tau0(l.tau0), H_MD(l.H_MD) {}

    // ! Destruction is automagic
    ~LatColMatPQPLeapfrogIntegrator(void) {};

    //! Get at the MD Hamiltonian
    AbsHamiltonian<multi1d<LatticeColorMatrix>,
		   multi1d<LatticeColorMatrix> >& getHamiltonian(void) {
      return *H_MD;
    }


  protected:

    void leapP(const Real& dt, 
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	                     multi1d<LatticeColorMatrix> >& s);
    //! Leap with Q
    void leapQ(const Real& dt, 
	       AbsFieldState<multi1d<LatticeColorMatrix>,
	                     multi1d<LatticeColorMatrix> >& s);

 

    //! Get the trajectory length
    const Real getTrajLength(void) const {
      return tau0;
    }

    //! Get the step size 
    const Real getStepSize(void) const {
      return dt;
    }

  private:
    Real dt;
    Real tau0;
    Handle< AbsHamiltonian<multi1d<LatticeColorMatrix>, 
			    multi1d<LatticeColorMatrix> > > H_MD;
    
    
  };

  





};


#endif
