#ifndef ILUPREC_S_CPREC_T_WILSONLIKE_LINOP_W_H
#define ILUPREC_S_CPREC_T_WILSONLIKE_LINOP_W_H

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
  
  class ILUPrecSCprecTWilsonLikeLinOp : 
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

    //! Destructor is automatic
    virtual ~ILUPrecSCprecTWilsonLikeLinOp() {}


  protected:
    

    virtual void spaceLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const = 0;
    virtual void derivSpaceOp(P& ds_u, const T& X, const T& Y, enum PlusMinus isign, int cb3d) const  = 0;

    virtual void AMinusOneOp(T& chi, const T& psi, enum PlusMinus isign) const =0 ;
    virtual void derivAMinusOne(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const = 0;

    //! Get log det ( T^\dag T )
    virtual Double logDetTDagT(void) const = 0;
    
    // Blanks to fill out...
    virtual const multi3d< int >&  getTsiteArray(void) const = 0;
    virtual const multi3d< CMat >& getPMatrixArray(void) const = 0;
    virtual const multi3d< CMat >& getPMatrixDaggerArray(void) const = 0;
    virtual const multi2d< CMat >& getQMatrixInvArray(void) const = 0;
    virtual const multi2d< CMat >& getQMatrixDaggerInvArray(void) const = 0;
    virtual const Real& getFactor() const = 0;
    virtual const Real& getInvFactor() const = 0;
    virtual const multi1d< LatticeColorMatrix>& getLinks(void) const = 0;
    virtual int getTMax(void) const = 0;
  public:

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const = 0;

    //! The time direction
    int tDir() const { 
      return 3; 
    }
    

    //! Apply inv (C_L)^{-1}
    virtual 
      void invCLeftLinOp(T& chi, const T& psi, enum PlusMinus isign) const
    { 

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

      getFermBC().modifyF(chi);

    }

    //! Apply inv (C_R)^{-1}
    virtual 
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
      getFermBC().modifyF(chi);
    }

    //! Apply C_L
    virtual 
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
      getFermBC().modifyF(chi);
    }
    
    //! Apply C_R
    virtual 
      void cRightLinOp(T& chi, const T& psi, enum PlusMinus isign) const {

      switch(isign) { 
      case PLUS: 
	{
	  //[ (T-)^{-1}_e                            0               ] 
	  //[ -(T-)^{-1}_{o} D_s_oe^{-1}   (T-)^{-1}_e (T-)^{-1}_o   ]      

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

      getFermBC().modifyF(chi);
    }

    // A = 0 so A-1 = -1 => (A-1) psi = -psi and also for the dagger.


    //! Flopcounter
    virtual unsigned long nFlops() const 
    { 
      //      QDPIO::cout << "Flopcount Not Yet Implemented " << endl;
      return 0;
    }

  protected:
    virtual
    void TPlusOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const
    {
      LatticeHalfFermion tmp_minus;
      LatticeHalfFermion tmp_plus;
      LatticeHalfFermion tmp_T;
      
      // This is just a way of doing P+ psi - decompose and reconstruct.
      // It may be more straightforward to just write a projector ... -- later
      tmp_plus[ rb3[ cb3d ] ]  = spinProjectDir3Plus(psi);
      chi[ rb3[ cb3d ] ] = spinReconstructDir3Plus(tmp_plus);

      // Here I project - apply TOp to only the halfector - then reco nstruct
      tmp_minus[ rb3[ cb3d ] ] = spinProjectDir3Minus(psi);


      // Use shared routine to apply T or T^\dagger.
      // Pass in u, and tsite for the subset
      CentralTPrecNoSpinUtils::TOp(tmp_T, 
				   tmp_minus, 
				   getLinks(),
				   getTsiteArray()[cb3d],
				   getFactor(),
				   isign,
				   schroedingerTP());

      chi[rb3[ cb3d ]]  += spinReconstructDir3Minus(tmp_T);
      chi[rb3[ cb3d ]]  *= Real(0.5);
      //      getFermBC().modifyF(chi);


    } 

    virtual
    void TMinusOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const 
    {
      enum PlusMinus other_sign = (isign == PLUS ? MINUS : PLUS) ;
      LatticeHalfFermion tmp_minus;
      LatticeHalfFermion tmp_plus;
      LatticeHalfFermion tmp_T;

      // This does the P- modulo a factor of 1/2
      // Rather than spooling through 2 half vectors I could just do a 
      // ProjectDir3Minus rather than going through half vectirs
      tmp_minus[rb3[cb3d]]  = spinProjectDir3Minus(psi);
      chi[rb3[cb3d]] = spinReconstructDir3Minus(tmp_minus);

      // This does the P+ T^\dagger modulo a factor of 1/2
      tmp_plus[rb3[cb3d]] = spinProjectDir3Plus(psi);


      // Use shared routine to apply T or T^\dagger.
      // Pass in u, and tsite for the subset
      CentralTPrecNoSpinUtils::TOp(tmp_T, 
				   tmp_plus, 
				   getLinks(),
				   getTsiteArray()[cb3d],
				   getFactor(),
				   other_sign,
				   schroedingerTP());

      chi[rb3[cb3d]] += spinReconstructDir3Plus(tmp_T);
      
      chi[rb3[cb3d]] *= Real(0.5); //The overall factor of 1/2

      //      getFermBC().modifyF(chi);


    }

    virtual
    void invTPlusOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const 
    {
      LatticeHalfFermion tmp_minus;
      LatticeHalfFermion tmp_plus;
      LatticeHalfFermion tmp_T;
      
      tmp_plus[rb3[cb3d]]  = spinProjectDir3Plus(psi);
      chi[rb3[cb3d]] = spinReconstructDir3Plus(tmp_plus);
      
      tmp_minus[rb3[cb3d]] = spinProjectDir3Minus(psi);
      
      // Use shared routine to apply (T)^{-1} or  (T^\dagger)^{-1}.
      // Pass in u, and tsite, P, P^\dag, Q and Q^\dag for the subset.
      CentralTPrecNoSpinUtils::invTOp( tmp_T,
				       tmp_minus, 
				       getLinks(),
				       getTsiteArray()[cb3d],
				       getPMatrixArray()[cb3d],
				       getPMatrixDaggerArray()[cb3d],
				       getQMatrixInvArray()[cb3d],
				       getQMatrixDaggerInvArray()[cb3d],
				       getInvFactor(),
				       isign,
				       getTMax(),
				       schroedingerTP());
      
      
      chi[rb3[cb3d]] += spinReconstructDir3Minus(tmp_T);
      chi[rb3[cb3d]] *= Real(0.5);
      //      getFermBC().modifyF(chi);


    }

    virtual
    void invTMinusOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const 
    {
      enum PlusMinus other_sign = isign == PLUS ? MINUS : PLUS ;
      LatticeHalfFermion tmp_minus;
      LatticeHalfFermion tmp_plus;
      LatticeHalfFermion tmp_T;
      
      tmp_minus[rb3[cb3d]]  = spinProjectDir3Minus(psi);
      chi[rb3[cb3d]] = spinReconstructDir3Minus(tmp_minus);
      
      tmp_plus[rb3[cb3d]] = spinProjectDir3Plus(psi);
      
      // Use shared routine to apply (T)^{-1} or  (T^\dagger)^{-1}.
      // Pass in u, and tsite, P, P^\dag, Q and Q^\dag for the subset.
      CentralTPrecNoSpinUtils::invTOp( tmp_T,
				       tmp_plus, 
				       getLinks(),
				       getTsiteArray()[cb3d],
				       getPMatrixArray()[cb3d],
				       getPMatrixDaggerArray()[cb3d],
				       getQMatrixInvArray()[cb3d],
				       getQMatrixDaggerInvArray()[cb3d],
				       getInvFactor(),	
				       other_sign,
				       getTMax(),
				       schroedingerTP());

      
      chi[rb3[cb3d]] += spinReconstructDir3Plus(tmp_T);
      chi[rb3[cb3d]] *= Real(0.5);
      //      getFermBC().modifyF(chi);


    }

    virtual
    void derivInvTPlusOp(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const 
    {
  
      // Same code as C_L in unprec case
      ds_u.resize(Nd);
      for(int mu=0; mu < 3; mu++) { 
	ds_u[mu] = zero;
      }
      
      if (isign == PLUS) { 
	
	LatticeHalfFermion tmp1, tmp2;
	T T1, T2; 
	
	tmp1=spinProjectDir3Minus(Y);
	for(int cb3d=0; cb3d < 2; cb3d++) { 
	  
	  CentralTPrecNoSpinUtils::invTOp( tmp2,
					   tmp1, 
					   getLinks(),
					   getTsiteArray()[cb3d],
					   getPMatrixArray()[cb3d],
					   getPMatrixDaggerArray()[cb3d],
					   getQMatrixInvArray()[cb3d],
					   getQMatrixDaggerInvArray()[cb3d],
					   getInvFactor(),	
					   PLUS,
					   getTMax(),
					   schroedingerTP());
	}
	
	T1  = spinReconstructDir3Minus(tmp2);
	
	tmp1=spinProjectDir3Minus(X);
	
	for(int cb3d=0; cb3d < 2; cb3d++) { 

	  CentralTPrecNoSpinUtils::invTOp( tmp2,
					   tmp1, 
					   getLinks(),
					   getTsiteArray()[cb3d],
					   getPMatrixArray()[cb3d],
					   getPMatrixDaggerArray()[cb3d],
					   getQMatrixInvArray()[cb3d],
					   getQMatrixDaggerInvArray()[cb3d],
					   getInvFactor(),
					   MINUS, 
					   getTMax(),
					   schroedingerTP());
	}
	
	// Two factors of 0.5 from the projectors.
	// Most cost efficient to apply them together to the half vector...
	tmp2 *= Real(0.25);
	
	T2  = spinReconstructDir3Minus(tmp2);
	
	
	LatticeFermion T3 = shift(T1, FORWARD, 3);
	
	// A minus from the fact that its dT^{-1}
	// A minus from the fact that the shift has a -
	// Makes  a + 
	ds_u[3] = traceSpin(outerProduct(T3, T2));
	
      }
      else {
	ds_u[3] = zero;
      }

      getFermBC().zero(ds_u);
    }

    virtual
    void derivInvTMinusOp(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const
    {
      // Same code as derivC_R 
      ds_u.resize(Nd);
      for(int mu=0; mu < 3; mu++) { 
	ds_u[mu] = zero;
      }
      
      if (isign == MINUS) { 
	
	LatticeHalfFermion tmp1, tmp2;
	T T1, T2; 
	
	tmp1=spinProjectDir3Plus(Y);
	for(int cb3d=0; cb3d < 2; cb3d++) { 
	  CentralTPrecNoSpinUtils::invTOp( tmp2,
					   tmp1, 
					   getLinks(),
					   getTsiteArray()[cb3d],
					   getPMatrixArray()[cb3d],
					   getPMatrixDaggerArray()[cb3d],
					   getQMatrixInvArray()[cb3d],
					   getQMatrixDaggerInvArray()[cb3d],
					   getInvFactor(),
					   PLUS,
					   getTMax(),
					   schroedingerTP());
	}      
	T1  = spinReconstructDir3Plus(tmp2);
	
	tmp1=spinProjectDir3Plus(X);
	for(int cb3d=0; cb3d < 2; cb3d++) { 
	  
	  CentralTPrecNoSpinUtils::invTOp( tmp2,
					   tmp1, 
					   getLinks(),
					   getTsiteArray()[cb3d],
					   getPMatrixArray()[cb3d],
					   getPMatrixDaggerArray()[cb3d],
					   getQMatrixInvArray()[cb3d],
					   getQMatrixDaggerInvArray()[cb3d],
					   getInvFactor(),
					   MINUS,
					   getTMax(),
					   schroedingerTP());
	}
	// Two factors of 0.5 from the projectors.
	// Most cost efficient to apply them together to the half vector...
	tmp2 *= Real(0.25);
	
	T2  = spinReconstructDir3Plus(tmp2);
	
	LatticeFermion T3 = shift(T1, FORWARD, 3);
	
	// A minus from the fact that its dT^{-1}
	// A minus from the fact that the shift has a -
	// Makes  a + 
	ds_u[3] = traceSpin(outerProduct(T3, T2));
	
      }
      else { 
	ds_u[3]= zero;
      }
      
      getFermBC().zero(ds_u);
    }
    
    virtual
      void  derivCRight(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const 
    {
      T T_1, T_2, T_3, T_4;
      switch(isign) {
      case PLUS: 
	{
	  // Simple Case: -T^o^\dagger_2 \dot(Dslash^oe_S) T^e_1
	  
	  // [ T^e_1 ]   [ (T-)^{-1}e  0 ] [ Y_e ]   [ (T-)^{-1}_e Y_e ]
	  // [       ] = [               ] [     ] = [                 ]
	  // { T^o_1 ]   [  0          1 ] [ Y_o ]   [   Y_o           ]
	  invTMinusOp(T_1, Y,  PLUS, 0);
	  T_1[ rb3[1] ] = Y;

	  // [ T^e_2 ]   [ 1           0   ] [ X_e ]  [   X_e             ]
	  // [       ] = [                 ] [     ] =[                   ]  
	  // [ T^o_2 ]   [ 0  (T-)^{-dag}o ] [ X_o ]  [ (T-)^{-dag}_o X_o ]
	  invTMinusOp(T_2, X, MINUS, 1);
	  T_2[ rb3[0] ] = X;
	  
	  derivSpaceOp(ds_u, T_2, T_1, PLUS, 1);
	  for(int mu=0; mu < Nd; mu++) { 
	    ds_u[mu] *= Real(-1);
	  }
	}
	break;
      case MINUS:
	{
	  // [ T^e_1 ]   [ 1           0   ] [ Y_e ]  [   Y_e             ]
	  // [       ] = [                 ] [     ] =[                   ]  
	  // [ T^o_1 ]   [ 0  (T-)^{-dag}o ] [ Y_o ]  [ (T-)^{-dag}_o Y_o ]
	  invTMinusOp(T_1, Y, MINUS, 1);
	  T_1[ rb3[0] ] = Y;
	  
	  
	  // [ T^e_2 ]   [ (T-)^{-1}e  0 ] [ X_e ]   [ (T-)^{-1}_e X_e ]
	  // [       ] = [               ] [     ] = [                 ]
	  // { T^o_2 ]   [  0          1 ] [ X_o ]   [   X_o           ]
	  invTMinusOp(T_2, X, PLUS, 0);
	  T_2[ rb3[1] ] = X;
	  
	  // [ T^e_3 ]   [ 1   - D^{oe}^\dag_s ][ T^e_1 ]
	  // [       ] = [                     ][       ]
	  // [ T^o_3 ]   [ 0          1        ][ T^o_1 ]
	  // We only need the even part tho, and put Y_o in the bottom
	  T T_tmp;
	  spaceLinOp(T_tmp, T_1, MINUS, 0);
	  T_3[rb3[0]] = T_1 - T_tmp; // Mhalf from the (-1/2)dslash
	  T_3[rb3[1]] = Y;
	  
	  
	  // [ T^e_4 ]   [ 1           0  ][ T^e_2 ]
	  // [       ] = [                ][       ]
	  // [ T^o_4 ]   [ -D^{eo}_s   1  ][ T^o_2 ]
	  // We only need the odd part tho, and put X_o in the top
	  spaceLinOp(T_tmp, T_2, PLUS, 1);
	  T_4[rb3[1]] = T_2 - T_tmp;       // Mhalf from the (-1/2)dslash
	  T_4[rb3[0]] = X;
	  
	  
	  derivInvTMinusOp(ds_u, T_4, T_3, MINUS);
	  
	  P ds_tmp;
	  
	  derivSpaceOp(ds_tmp, T_2, T_1, MINUS, 0);
	  for(int mu=0; mu < Nd; mu++) {
	    ds_u[mu] -= ds_tmp[mu];           
	  }
	  
	}
	break;
      default:
	QDPIO::cerr << "Bad Case. Should never get here" << endl;
	QDP_abort(1);
      }
      getFermBC().zero(ds_u);
    }

    virtual 
    void  derivCLeft(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const
    {
      T T_1, T_2, T_3, T_4;
      
      switch(isign) {
      case PLUS: 
	{
	  // [ T^e_1 ]   [ 1      0        ] [ Y_e ]  [   Y_e             ]
	  // [       ] = [                 ] [     ] =[                   ]  
	  // [ T^o_1 ]   [ 0   (T+)^{-1}_o ] [ Y_o ]  [ (T+)^{-1}_o Y_o   ]
	  invTPlusOp(T_1, Y, PLUS, 1);
	  T_1[ rb3[0] ] = Y;
	  
	  
	  // [ T^e_2 ]   [ (T+)^{-dag}e  0 ] [ X_e ]   [ (T+)^{-dag}_e X_e ]
	  // [       ] = [                 ] [     ] = [                   ]
	  // { T^o_2 ]   [  0            1 ] [ X_o ]   [   X_o             ]
	  invTPlusOp(T_2, X, MINUS, 0);
	  T_2[ rb3[1] ] = X;
	  
	  // [ T^e_3 ]   [ 1   - D^{eo}_s ][ T^e_1 ]
	  // [       ] = [                ][       ]
	  // [ T^o_3 ]   [ 0          1   ][ T^o_1 ]
	  // We only need the even part tho, and put Y_o in the bottom
	  T T_tmp;
	  spaceLinOp(T_tmp, T_1, PLUS, 0);
	  T_3[rb3[0]] = T_1 - T_tmp;  // mhalf is from (-1/2)Dslash
	  T_3[rb3[1]] = Y;
	  
	  
	  // [ T^e_4 ]   [ 1                0  ][ T^e_2 ]
	  // [       ] = [                     ][       ]
	  // [ T^o_4 ]   [ -D^{oe \dag}_s   1  ][ T^o_2 ]
	  // We only need the odd part tho, and put X_o in the top
	  spaceLinOp(T_tmp, T_2, MINUS, 1);
	  T_4[rb3[1]] = T_2 - T_tmp; // mhalf is from (-1/2) Dslash
	  T_4[rb3[0]] = X;
	  
	  P ds_tmp;
	  derivSpaceOp(ds_tmp, T_2, T_1, PLUS, 0);
	  derivInvTPlusOp(ds_u, T_4, T_3, PLUS);
	  for(int mu=0; mu < Nd; mu++) { 
	    ds_u[mu] -= ds_tmp[mu];
	  }
	}
	break;
	
      case MINUS:
	{
	  // Simple Case: -T^o^\dagger_2 \dot(Dslash^oe_S) T^e_1
	  
	  // [ T^e_1 ]   [ (T+)^{-\dag}e  0 ] [ Y_e ]   [ (T+)^{-\dag}_e Y_e ]
	  // [       ] = [                  ] [     ] = [                 ]
	  // { T^o_1 ]   [  0             1 ] [ Y_o ]   [   Y_o           ]
	  invTPlusOp(T_1, Y,  MINUS, 0);
	  T_1[ rb3[1] ] = Y;
	  
	  // [ T^e_2 ]   [ 1           0   ] [ X_e ]  [   X_e             ]
	  // [       ] = [                 ] [     ] =[                   ]  
	  // [ T^o_2 ]   [ 0  (T+)^{-1}o   ] [ X_o ]  [ (T+)^{-1}_o X_o ]
	  invTPlusOp(T_2, X, PLUS, 1);
	  T_2[ rb3[0] ] = X;
	  
	  derivSpaceOp(ds_u, T_2, T_1, MINUS, 1);
	  for(int mu=0; mu < Nd; mu++) { 
	    ds_u[mu] *= Real(-1);  // (Combinme factor of -1 and -(1/2)
	  }
	}
	break;
      default:
	QDPIO::cerr << "Bad Case. Should never get here" << endl;
	QDP_abort(1);
      }

      getFermBC().zero(ds_u);      
    }
    
    //! Get the force due to the det T^\dag T bit
    //  This is the same as for the unprec case as it turns
    virtual void derivLogDetTDagT(P& ds_u, enum PlusMinus isign) const 
    {
      // Derivative of a Hermitian quantity so ignore isign?
      // Initial development -- set to 0
      ds_u.resize(Nd);
	

      for(int mu=0; mu < Nd; mu++) { 
	ds_u[mu] = zero;
      }

      if( !schroedingerTP() ) {

	// Force from even sites
	CentralTPrecNoSpinUtils::derivLogDet(ds_u, 
					     getLinks(),
					     getQMatrixInvArray(),
					     getTsiteArray(),
					     tDir(),
					     getFactor(),
					     schroedingerTP());
	
	
	
	getFermBC().zero(ds_u);
      }
      else {
	ds_u[tDir()] = zero;
      }
    }

  };



} // End Namespace Chroma

#endif

#endif
#endif
#endif
