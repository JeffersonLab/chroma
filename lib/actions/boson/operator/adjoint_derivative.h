// -*- C++ -*-
/*! \file
 *  \brief adjoint derivative 
 */

#ifndef __adjoint_derivative_h__
#define __adjoint_derivative_h__

#include "linearop.h"
#include "util/gauge/stag_phases_s.h"

namespace Chroma 
{ 
  //! Adjoint derivative operator
  /*!
   * \ingroup linop
   *
   * ADJOINT_DERIVATIVE
   *
   */
  class AdjointDerivative : 
    public LinearOperator<LatticeColorMatrix>
  {
  public:
    AdjointDerivative(int mu_, Real rho_, multi1d<LatticeColorMatrix>& u_):mu(mu_),rho(rho_), u(u_) {
      // Assume periodic gauge BC
      // add the staggered phases...
      for(int i = 0; i < Nd; i++) 
	{
	  u[i]     *= StagPhases::alpha(i);
	}
      halfI=Complex(0.0,0.5); // define i/2 for optimization purposes
    }
    //! No real need for cleanup here
    ~AdjointDerivative() {}

    //! Subset is all here
    const Subset& subset() const {return all;}

    //! Return flops performed by the operator()
    // approximatelly correct
    unsigned long nFlops(){return 4*Nc*Nc + (Nd-1)*(4*Nc*Nc*Nc + 2Nc*Nc) ;  } const;

    virtual void operator()(LatticeColorMatrix& chi, const LatticeColorMatrix psi, enum PlusMinus ising) const{
      // Implement the derivative
      //this is a hermitian operator therefore PlusMinus makes no difference
      LatticeColorMatrix foo = zero ;
      for(int nu(0);nu  < Nd; nu++ ){
	if(nu!=mu){
	  foo = foo +
	    u[nu]*shift(psi, FORWARD, nu)*adj(u[nu]) -
	    shift(adj(u[nu])*psi*u[nu],BACKWARD,nu) ;
	}
      }
      chi = rho*psi + halfI*foo ;
    }
  protected:

  private:
    Real rho ; 
    int mu ; // if mu in [0, Nd-1] do the transverse Nd-1 gradient
    multi1d<LatticeColorMatrix> u ; // u with phases
    Complex halfI ;
  };

  class squaredAdjointDerivative : 
    public LinearOperator<LatticeColorMatrix>
  {
  public:
    squaredAdjointDerivative(const AdjointDerivative& D_):D(D_) {}
    //! No real need for cleanup here
    ~squaredAdjointDerivative() {}
 //! Subset is all here
    const Subset& subset() const {return all;}

    //! Return flops performed by the operator()
    unsigned long nFlops(){return 2*D.hFlops } const;

    virtual void operator()(LatticeColorMatrix& chi, const LatticeColorMatrix psi, enum PlusMinus ising) const{
      // Implement the derivative
      //this is a hermitian operator therefore PlusMinus makes no difference
      LatticeColorMatrix foo  ;
      D(foo,psi);
      D(chi,foo);
    }
  protected:

  private:
    AdjointDerivative D ;
    
  };
} // End Namespace Chroma


#endif
