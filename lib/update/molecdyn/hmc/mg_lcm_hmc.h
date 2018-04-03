// -*- C++ -*-
/*! \file
 * \brief HMC trajectory
 *
 * HMC trajectory
 */

// Implements the hyperplane multi-scale HMC

#ifndef MG_LAT_COL_MAT_HMC_NEW_H
#define MG_LAT_COL_MAT_HMC_NEW_H

#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/hamiltonian/abs_hamiltonian.h"
#include "update/molecdyn/integrator/abs_integrator.h"
#include "update/molecdyn/hmc/abs_hmc.h"
#include "handle.h"
// The accept Reject
#include "update/molecdyn/hmc/global_metropolis_accrej.h"

#include "util/gauge/taproj.h" 

#include "actions/boson/operator/adjoint_derivative.h"

namespace Chroma 
{

  //! HMC trajectory
  /*! @ingroup hmc */
  class MGLatColMatHMCTrj : public AbsMGHMCTrj<multi1d<LatticeColorMatrix>, 
			  multi1d<LatticeColorMatrix> > 
  {
  public:

    // Destructor
    ~MGLatColMatHMCTrj(void) {} 

    // Constructor
    MGLatColMatHMCTrj( Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& _H_MC,
			Handle< AbsMDIntegrator< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > >& _MD_int)
      : the_MD(_MD_int), the_H_MC(_H_MC),dir(-1),rho(0.01) {}

  private:
    Handle< AbsMDIntegrator<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > the_MD; 

    Handle< AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > the_H_MC;

    int dir ;
    Real rho ;
    
  protected:

    AbsHamiltonian< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& getMCHamiltonian(void) { 
      return *the_H_MC;
    }

    AbsMDIntegrator< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& getMDIntegrator(void) {
      return *the_MD;
    }

    bool acceptReject( const Double& DeltaH ) const {
      return globalMetropolisAcceptReject(DeltaH);
    }

    void setDir(int d){dir=d;}
    const int& theDir() const {return dir;}

    void setRho(Real d){rho=d;}
    const Real& theRho() const {return rho;} 
    
    void refreshP(AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
    {
      START_CODE();

      
      if((dir>-1)&&(dir<Nd)){//doing new algorithm
	// Loop over direcsions
	
	for(int mu = 0; mu < Nd; mu++) 
	  {
	    s.getP()[mu] = zero ;
	    if(mu==dir){
	      LatticeColorMatrix P ;
	      gaussian(P);
	      P *= sqrt(Real(0.5));
	      // Make traceless and antihermitian
	      taproj(P);
	      
	      AdjointDerivative D(mu,rho,s.getQ());
	      D(s.getP()[mu],P,MINUS) ;
	      // I need to check here that inded the momenta are traceless
	      // anti hermitian
	    }
	
	  }
      }
      else{// do normal HMC
	// Loop over direcsions
	for(int mu = 0; mu < Nd; mu++) 
	  {
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
    
      END_CODE();
    }


    void flipMomenta(AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const {

      multi1d<LatticeColorMatrix>& p = s.getP();
      for(int mu=0; mu < Nd; mu++) { 
	p[mu] *= Real(-1);
      }
    }

    void reverseCheckMetrics(Double& deltaQ, Double& deltaP,
			     const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s,
			     const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s_old
			     ) const  { 

      multi1d<LatticeColorMatrix> DelQ(Nd);
      multi1d<LatticeColorMatrix> DelP(Nd);
      for(int mu=0; mu < Nd; mu++) {
	DelQ[mu]=s.getQ()[mu] - s_old.getQ()[mu];
	DelP[mu]=s.getP()[mu] - s_old.getP()[mu];
      }
      
      deltaQ = sqrt(norm2(DelQ[0]));
      deltaP = sqrt(norm2(DelP[0]));
      for(int mu=1; mu < Nd; mu++) { 
	deltaQ += sqrt(norm2(DelQ[mu]));
	deltaP += sqrt(norm2(DelP[mu]));
      }

      deltaQ /= Double(Nd*Layout::vol());
      deltaP /= Double(Nd*Layout::vol());

      
    }

  };
}

#endif
