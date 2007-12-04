#ifndef ILUPREC_S_CPREC_T_CLOVER_LINOP_W_H
#define ILUPREC_S_CPREC_T_CLOVER_LINOP_W_H
#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3

#include "linearop.h"
#include "central_tprec_linop.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/central_tprec_nospin_utils.h"
#include "actions/ferm/linop/clover_term_w.h"

namespace Chroma 
{ 
  //! Clover Dirac Operator - Unpreconditioned in Space, Centrally Preconditioned in time
  /*!
   * \ingroup linop
   *
   * This routine is specific to Clover fermions!
   *
   *                                                      ~      ~+
   */
  
  class ILUPrecSCprecTCloverLinOp : 
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
    ILUPrecSCprecTCloverLinOp() {}


    //! Full constructor with Anisotropy
    ILUPrecSCprecTCloverLinOp(Handle< FermState<T,P,Q> > fs_,
			      const CloverFermActParams& param_ )
      {create(fs_,param_);}

    //! Destructor is automatic
    ~ILUPrecSCprecTCloverLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return fs->getBC();}


    //! Creation routine with Anisotropy
    void create(Handle< FermState<T,P,Q> > fs_,
		const CloverFermActParams& param_);


    //! The time direction
    
    int tDir() const { 
      // Always 3 for now?
      return 3; 
    }
    
    void deriv(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const {
      QDPIO::cout << "Not implemented " << endl;

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

    //! Do chi = (A - 1) psi with A = -(c_sw/4) sigma_munu F_munu
    inline
    void AMinusOneOp(T& chi, const T& psi, enum PlusMinus isign) const { 
      T tmp1;

      // The clover term does A_plus_fact = factI - (c_sw/4) sigma_munu F_munu = A + fact
      // where the desired A is just A = - (c_sw/4) sigma_munu F_munu

      // I need to get rid of the factI. Luckily, since I already
      // have to do A - 1, I can just change this to A-(1 + fact)I

      // 1 + fact
      Real ftmp = Real(1) + fact;

      // tmp = (A + fact) psi
      APlusFact(tmp1, psi, isign);

      // chi
      chi = tmp1 - ftmp*psi;
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


  private:


    Real fact;  // tmp holding  Nd+Mass
    Real invfact;
    Handle< FermState<T,P,Q> > fs;
    multi1d<LatticeColorMatrix> u;

    multi3d< int  > tsite; // Tsite array - 3d even sites
       
    multi3d< CMat > P_mat; // 3d even sites
    multi3d< CMat > P_mat_dag;

    multi2d< CMat > Q_mat_inv;        // Just one matrix for each spatial site ( 1+P[t=0] )^{-1}
    multi2d< CMat > Q_mat_dag_inv;    // Just one matrix for each spatial site ( 1+P_dag[t=Nt-1])^{-1}


    CloverTerm          APlusFact;
    CloverFermActParams param;

    WilsonDslash3D Dw3D;
  };

} // End Namespace Chroma

#endif
#endif
#endif

#endif
