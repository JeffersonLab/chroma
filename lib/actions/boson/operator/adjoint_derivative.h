// -*- C++ -*-
/*! \file
 *  \brief adjoint derivative 
 */

#ifndef __adjoint_derivative_h__
#define __adjoint_derivative_h__

#include "linearop.h"

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
    // I need to check if the staggered phases work if restricted to 3D
  public:
    AdjointDerivative(int mu_, Real rho_, const multi1d<LatticeColorMatrix>& u_):mu(mu_),rho(rho_), u(u_) {
      if(Nd!=4){
	QDPIO::cout<<"AdjointDerivative: works only for Nd=4. We have Nd ="
		   <<Nd<<std::endl ;
	QDP_abort(2923);
      }
      // Assume periodic gauge BC
      // add the staggered phases...
      multi1d<LatticeInteger> phases(Nd);

      multi1d<LatticeInteger> x(Nd);
      // Fill x with lattice coordinates
      for( int nu = 0; nu < Nd; nu++) {
        x[ nu ] = Layout::latticeCoordinate(nu);
      }
      
      phases[0] = LatticeInteger(1);
      
    
      switch (mu){
      case 0:
	phases[1] = LatticeInteger(1);
	phases[2]=where( x[1] % 2 == 0, LatticeInteger(1), LatticeInteger(-1));
	phases[3]=where( ((x[1]+x[2])%2) == 0, LatticeInteger(1), LatticeInteger(-1));
	break ;
      case 1:
	phases[1] = LatticeInteger(1);
	phases[2]=where( x[0] % 2 == 0, LatticeInteger(1), LatticeInteger(-1));
	phases[3]=where( ((x[0]+x[2])%2) == 0, LatticeInteger(1), LatticeInteger(-1));
	break ;
      case 2:
	phases[1]=where( x[0] % 2 == 0, LatticeInteger(1), LatticeInteger(-1));
	phases[2] = LatticeInteger(1);
	phases[3]=where( ((x[0]+x[1])%2) == 0, LatticeInteger(1), LatticeInteger(-1));
	break ;
      case 3:
	phases[1]=where( x[0] % 2 == 0, LatticeInteger(1), LatticeInteger(-1));
	phases[2]=where( ((x[0]+x[1])%2) == 0, LatticeInteger(1), LatticeInteger(-1));
	phases[3]=LatticeInteger(1);
	break ;
      default:
	// do the 4-d case
	//phases[0] = LatticeInteger(1);
        phases[1] = where( x[0] % 2 == 0, LatticeInteger(1), LatticeInteger(-1));
        phases[2] = where( (x[0]+x[1] ) % 2  == 0, LatticeInteger(1), LatticeInteger(-1));
        phases[3] = where( (x[0]+x[1]+x[2] ) % 2  == 0, LatticeInteger(1), LatticeInteger(-1));
	break ;
      }

      for(int i = 0; i < Nd; i++) 
	{
	  u[i]     *= phases[i];
	}
      halfI=cmplx(Real(0.0),Real(0.5)); // define i/2 for optimization purposes
    }
    //! No real need for cleanup here
    ~AdjointDerivative() {}

    //! Subset is all here
    const Subset& subset() const {return all;}

    //! Return flops performed by the operator()
    // approximatelly correct
    unsigned long nFlops() const {
      return 4*Nc*Nc + (Nd-1)*(4*Nc*Nc*Nc + 2*Nc*Nc) ;
    }

    virtual void operator()(LatticeColorMatrix& chi, const LatticeColorMatrix& psi, enum PlusMinus ising) const {
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

      // the derivative is anti-hermitian
      // isign negative means the dagger 
      Real half = 0.5*ising ;
      chi = rho*psi + half*foo ;
	  
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
    unsigned long nFlops() const {return 2*D.nFlops(); } 

    virtual void operator()(LatticeColorMatrix& chi, const LatticeColorMatrix& psi, enum PlusMinus ising) const {
      // Implement the derivative
      //this is a hermitian operator therefore PlusMinus makes no difference
      LatticeColorMatrix foo  ;
      D(foo,psi,PLUS);
      D(chi,foo,MINUS);
    }
  protected:

  private:
    AdjointDerivative D ;
    
  };
} // End Namespace Chroma


#endif
