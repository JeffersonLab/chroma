#ifndef ILUPREC_S_CPREC_T_WILSON_LINOP_W_H
#define ILUPREC_S_CPREC_T_WILSON_LINOP_W_H

#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3

#include "linearop.h"
#include "central_tprec_linop.h"

#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/central_tprec_nospin_utils.h"

namespace Chroma 
{ 
  //! Wilson Dirac Operator - Unpreconditioned in Space, Centrally Preconditioned in time
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   *                                                      ~      ~+
   */
  
  class ILUPrecSCprecTWilsonLinOp : 
    public ILUPrecSpaceCentralPrecTimeLinearOperator<LatticeFermion, 
    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;
    typedef PScalar<PColorMatrix<RComplex<REAL>, 3> > CMat;              // Useful type: ColorMat with no Outer<>
    typedef PSpinVector<PColorVector<RComplex<REAL>, 3>, 2> HVec_site;   // Useful type: Half Vec with no Outer<>
    //! Partial constructor
    ILUPrecSCprecTWilsonLinOp() {}


    //! Full constructor with Anisotropy
    ILUPrecSCprecTWilsonLinOp(Handle< FermState<T,P,Q> > fs_,
		      const Real& Mass_,
		      const AnisoParam_t& aniso_)
      {create(fs_,Mass_,aniso_);}

    //! Destructor is automatic
    ~ILUPrecSCprecTWilsonLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return fs->getBC();}


    //! Creation routine with Anisotropy
    void create(Handle< FermState<T,P,Q> > fs_,
		const Real& Mass_,
		const AnisoParam_t& aniso_);


    //! The time direction
    
    int tDir() const { 
      // Always 3 for now?
      return 3; 
    }
    

    //! Apply inv (C_L)^{-1}
    inline
    void invCLeftLinOp(T& chi, const T& psi, enum PlusMinus isign) const { 

      switch(isign) { 
      case PLUS: 
	{
	  //[ T+   D_s ] = 
	  //[ 0     T+ ] 
	  T tmp1;
	  TPlusOp(chi, psi, PLUS, 0);
	  TPlusOp(chi, psi, PLUS, 1);
	  spaceLinOp(tmp1, psi, PLUS, 0);
	  chi[rb3[0]] += tmp1;
	}
	break;
      case MINUS:
	{
	  //
	  // Remember to dagger in Even Odd Space too.
	  // [  (T+)^\dag   0         ] 
	  // [    D_s^\dag  (T+)^\dag ]
	  T tmp1;
	  TPlusOp(chi, psi, MINUS, 0);
	  TPlusOp(chi, psi, MINUS, 1);
	  spaceLinOp(tmp1, psi, MINUS, 1);
	  chi[rb3[1]] += tmp1;
	}
	break;
      default:
	QDPIO::cout << "Unknown sign " << endl;
	QDP_abort(1);
      }
    }

    //! Apply inv (C_R)^{-1}
    inline
    void invCRightLinOp(T& chi, const T& psi, enum PlusMinus isign) const {
      switch(isign) { 
      case PLUS: 
	{
	  //[ T-      0  ] 
	  //[ D_s     T- ]
	  T tmp1;
	  TMinusOp(chi, psi, PLUS, 0);
	  TMinusOp(chi, psi, PLUS, 1);
	  spaceLinOp(tmp1, psi, PLUS, 1);
	  chi[rb3[1]] += tmp1;
	}
	break;
      case MINUS:
	{
	  //
	  // Remember to dagger in Even Odd Space too.
	  // [  (T-)^\dag     D_s^\dag ] 
	  // [       0       (T-)^\dag  ]
	  T tmp1;
	  TMinusOp(chi, psi, MINUS, 0);
	  TMinusOp(chi, psi, MINUS, 1);
	  spaceLinOp(tmp1, psi, MINUS, 0);
	  chi[rb3[0]] += tmp1;
	}
	break;
      default:
	QDPIO::cout << "Unknown sign " << endl;
	QDP_abort(1);
      };
    }

    //! Apply C_L
    inline
    void cLeftLinOp(T& chi, const T& psi, enum PlusMinus isign) const { 
      switch(isign) { 
      case PLUS: 
	{
	  //[ (T+)^{-1}  -(T+)^{-1}_{e} D_s_eo (T+)^{-1}_o ] 
	  //[ 0                          (T+)^{-1}_o       ] 
	  T tmp1, tmp2;
	  

	  // tmp1_0 = (T+)^{-1}_{o} \psi_o
	  invTPlusOp(tmp1, psi, PLUS, 1);

	  // chi_o = (T+)^{-1}_{o} \psi_o
	  chi[rb3[1]] = tmp1;

	  // tmp2_e = D_s_eo (T+)^{-1}_{o} \psi_o
	  spaceLinOp(tmp2, tmp1, PLUS, 0);

	  // tmp1_e = \psi_e - D_s_eo (T+)^{-1}_{o} \psi_o
	  tmp1[rb3[0]] = psi - tmp2;


	  // chi_e = (T+)^{-1}_{e} ( \psi_e -  D_s_eo (T+)^{-1}_{o} \psi_o ) 
	  invTPlusOp(chi, tmp1, PLUS, 0);

	}
	break;
      case MINUS:
	{
	  //[ (T+)^{-\dag                                   0               ] 
	  //[ -(T+)^{-\dag}_{o} D_s_oe^\dag (T+)^{-\dag}_e (T+)^{-\dag}_o   ]      

	  // = [ 1     0          ] [ T+^{-\dag}              0 ] 
	  //   [ 0   T+_o^{-\dag} ] [ -D_{oe}^\dag T+^{-\dag}  1 ]
	  T tmp1, tmp2;
	  
	  // tmp1 = T+^{-dag}_{e} \psi_{e}
	  invTPlusOp(tmp1, psi, MINUS, 0);

	  // chi_e =  T+^{-dag}_{e} \psi_{e}
	  chi[rb3[0]] = tmp1;

	  // tmp2_o = D_{oe} T+^{-\dag} psi_e 
	  spaceLinOp(tmp2, tmp1, MINUS, 1);

	  // tmp1_o = psi_o -  D_{oe} T+^{-\dag} psi_e 
	  tmp1[rb3[1]] = psi - tmp2;

	  invTPlusOp(chi, tmp1, MINUS, 1);
	}
	break;
      default:
	QDPIO::cout << "Unknown sign " << endl;
	QDP_abort(1);
      };
    }
    
    //! Apply C_R
    inline
    void cRightLinOp(T& chi, const T& psi, enum PlusMinus isign) const {
      switch(isign) { 
      case PLUS: 
	{
	  //[ (T-)^{-1}_e                            0               ] 
	  //[ -(T-)^{-1}_{o} D_s_oe^{-1} (T-)^{-1}_e (T-)^{-1}_o   ]      

	  // = [ 1     0       ] [ T-^{-1}_e              0 ] 
	  //   [ 0   T-_o^{-1} ] [ -D_{oe}^\dag T-^{-1}_e  1 ]
	  T tmp1, tmp2;
	  
	  // tmp1 = T+^{-dag}_{e} \psi_{e}
	  invTMinusOp(tmp1, psi, PLUS, 0);

	  // chi_e =  T+^{-dag}_{e} \psi_{e}
	  chi[rb3[0]] = tmp1;

	  // tmp2_o = D_{oe} T+^{-\dag} psi_e 
	  spaceLinOp(tmp2, tmp1, PLUS, 1);

	  // tmp1_o = psi_o -  D_{oe} T+^{-\dag} psi_e 
	  tmp1[rb3[1]] = psi - tmp2;

	  invTMinusOp(chi, tmp1, PLUS, 1);
	}


	break;
      case MINUS:
      	{
	  //[ (T-)^{-dag}  -(T-)^{-\dag}_{e} D_s_eo^\dag (T-)^{-\dag}_o ] 
	  //[ 0                          (T-)^{-\dag}_o                 ] 
	  T tmp1, tmp2;
	  

	  // tmp1_0 = (T-)^{-\dag}_{o} \psi_o
	  invTMinusOp(tmp1, psi, MINUS, 1);

	  // chi_o = (T-)^{-\dag}_{o} \psi_o
	  chi[rb3[1]] = tmp1;

	  // tmp2_e = D_s_eo (T-)^{-1\dag_{o} \psi_o
	  spaceLinOp(tmp2, tmp1, MINUS, 0);

	  // tmp1_e = \psi_e - D_s_eo (T-)^{-\dag}_{o} \psi_o
	  tmp1[rb3[0]] = psi - tmp2;


	  // chi_e = (T-)^{-\dag}_{e} ( \psi_e -  D_s_eo (T-)^{-\dag}_{o} \psi_o ) 
	  invTMinusOp(chi, tmp1, MINUS, 0);

      }
	break;
      default:
	QDPIO::cout << "Unknown sign " << endl;
	QDP_abort(1);
      };
    }

    // A = 0 so A-1 = -1 => (A-1) psi = -psi and also for the dagger.
    inline
    void AMinusOneOp(T& chi, const T& psi, enum PlusMinus isign) const { 
      chi = Real(-1)*psi;
    }
    
    //! Flopcounter
    unsigned long nFlops() const 
    { 
      //      QDPIO::cout << "Flopcount Not Yet Implemented " << endl;
      return 0;
    }
    
  private:
    
    //! Apply the the space block onto a source vector
    //  cb3d is the 3d (rb3) checkerboard of the target
    void spaceLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const {
      Real mhalf=Real(-0.5);
      Dw3D.apply(chi, psi, isign, cb3d);
      chi[ rb3[cb3d] ] *= mhalf;
    }

    void TPlusOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const;

    void TMinusOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const;

    void invTPlusOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const;

    void invTMinusOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const;



    //! Apply the d/dt of the preconditioned linop
    void deriv(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const {

      QDPIO::cerr << "Not Yet Implemented" << endl;
      QDP_abort(1);

    }
    
    //! Get log det ( T^\dag T )
    Double logDetTDagT(void) const {
      return logDetTSq;
    }


    //! Get the force due to the det T^\dag T bit
    void derivLogDetTDagT(P& ds_u, enum PlusMinus isign) const {
      QDPIO::cerr << "Not Yet Implemented" << endl;
      QDP_abort(1);
    }


  private:

    AnisoParam_t aniso;
    Real fact;  // tmp holding  Nd+Mass
    Real invfact;
    Real Mass;
    Handle< FermState<T,P,Q> > fs;
    multi1d<LatticeColorMatrix> u;

    multi3d< int  > tsite; // Tsite array - 3d even sites
       
    multi3d< CMat > P_mat; // 3d even sites
    multi3d< CMat > P_mat_dag;

    multi2d< CMat > Q_mat_inv;        // Just one matrix for each spatial site ( 1+P[t=0] )^{-1}
    multi2d< CMat > Q_mat_dag_inv;    // Just one matrix for each spatial site ( 1+P_dag[t=Nt-1])^{-1}
    
    WilsonDslash3D Dw3D;
    Double logDetTSq;
  };

} // End Namespace Chroma

#endif

#endif
#endif
#endif
