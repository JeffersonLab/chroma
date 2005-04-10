// -*- C++ -*-
// $Id: lcm_hmc.h,v 1.4 2005-04-10 21:46:42 edwards Exp $
/*! \file
 * \brief HMC trajectory
 *
 * HMC trajectory
 */

#ifndef LAT_COL_MAT_HMC_H
#define LAT_COL_MAT_HMC_H

#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/hmc/abs_hmc.h"

// The accept Reject
#include "update/molecdyn/hmc/global_metropolis_accrej.h"

#include "util/gauge/taproj.h" 


namespace Chroma 
{

  //! HMC trajectory
  /*! @ingroup hmc */
  class LatColMatHMCTrj : public AbsHMCTrj<multi1d<LatticeColorMatrix>, 
			  multi1d<LatticeColorMatrix> > {
  public:

    // Destructor
    ~LatColMatHMCTrj(void) {} 

    // Constructor
    LatColMatHMCTrj( Handle< ExactAbsHamiltonian< multi1d<LatticeColorMatrix>, 
		     multi1d<LatticeColorMatrix> > >& _H_MC,
		     Handle< AbsMDIntegrator< multi1d<LatticeColorMatrix>, 
		     multi1d<LatticeColorMatrix> > >& _MD_integrator ) : the_MD(_MD_integrator), the_H_MC(_H_MC) {}

  private:
    Handle< AbsMDIntegrator<multi1d<LatticeColorMatrix>, 
			    multi1d<LatticeColorMatrix> > > the_MD;

    Handle< ExactAbsHamiltonian< multi1d<LatticeColorMatrix>, 
				 multi1d<LatticeColorMatrix> > > the_H_MC;

  protected:

    ExactAbsHamiltonian< multi1d<LatticeColorMatrix>, 
			 multi1d<LatticeColorMatrix> >& getMCHamiltonian(void) { 
      return *the_H_MC;
    }

    AbsMDIntegrator< multi1d<LatticeColorMatrix>, 
		     multi1d<LatticeColorMatrix> >& getMDIntegrator(void) {
      return *the_MD;
    }

    bool acceptReject( const Double& DeltaH ) const {
      return globalMetropolisAcceptReject(DeltaH);
    }


    void refreshP(AbsFieldState<multi1d<LatticeColorMatrix>, 
		  multi1d<LatticeColorMatrix> >& s) const {
      
      // Loop over direcsions
      for(int mu = 0; mu < Nd; mu++) {

	// Pull the gaussian noise
	gaussian(s.getP()[mu]);

	// Old conventions
	//s.getP()[mu] *= sqrt(0.5);  // Gaussian Normalisation

	// Normalisation factor in variance of momenta 
	// one factor of sqrt(2) to move the variance of the noise
	// one factor of sqrt(2) to account for Taproj normalisation
	// => sqrt(1/2)*sqrt(1/2) = sqrt(1/4) = 1/2 
	s.getP()[mu] *= sqrt(Real(0.5)); 
                         
	// Make traceless and antihermitian
	taproj(s.getP()[mu]);

      }
    }
  };
};

#endif
