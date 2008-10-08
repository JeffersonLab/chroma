#ifndef ILU2PREC_S_CPREC_T_WILSONLIKE_LINOP_W_H
#define ILU2PREC_S_CPREC_T_WILSONLIKE_LINOP_W_H

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
  
  class ILU2PrecSCprecTWilsonLikeLinOp : 
    public ILU2PrecSpaceCentralPrecTimeLinearOperator<LatticeFermion, 
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
    virtual ~ILU2PrecSCprecTWilsonLikeLinOp() {}


  protected:
    

    
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
    

    // A = 0 so A-1 = -1 => (A-1) psi = -psi and also for the dagger.


    //! Flopcounter
    virtual unsigned long nFlops() const 
    { 
      //      QDPIO::cout << "Flopcount Not Yet Implemented " << endl;
      return 0;
    }

  protected:

    // C_L^{-1} P+ + P_ T
    virtual
    void invCLeftLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const
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


    // C^{-1}_R =  P- + P+ T^\dagger
    virtual
    void invCRightLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const 
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


    // C_L = P+ + P_ T^{-1}
    virtual
    void cLeftLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const 
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


    // C_R =  P- + P+ T^{-\dagger}
    virtual
    void cRightLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const 
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

    virtual void DBar(T& chi, const T& psi, enum PlusMinus isign, int cb3) const
    {
      T tmp1, tmp2;
      if (isign == PLUS) { 
	cRightLinOp(tmp1, psi, PLUS, (1-cb3));
	Dslash3D(tmp2, tmp1, PLUS, cb3);
	cLeftLinOp(chi, tmp2, PLUS, cb3);
      }
      else {
	cLeftLinOp(tmp1, psi, MINUS, 1-cb3);
	Dslash3D(tmp2, tmp1, MINUS, cb3);
	cRightLinOp(chi, tmp2, MINUS, cb3);
      }
    }

     
    virtual void ABar(T& chi, const T& psi, enum PlusMinus isign, int cb3) const
    {
      T tmp1, tmp2; 
      if( isign == PLUS ) { 
	cRightLinOp(tmp1, psi, PLUS, cb3);
	AH(tmp2, tmp1, PLUS, cb3);
	cLeftLinOp(chi, tmp2, PLUS, cb3);
      }
      else {
	cLeftLinOp(tmp1, psi, MINUS, cb3);
	AH(tmp2, tmp1, MINUS, cb3);
	cRightLinOp(chi, tmp2, MINUS, cb3);
      }
    }

    virtual void Dslash3D(T& chi, const T& psi, enum PlusMinus isign, int cb) const = 0;
    virtual void AH(T& chi, const T& psi, enum PlusMinus isign, int cb) const = 0;

  };


} // End Namespace Chroma

#endif

#endif
#endif
#endif
