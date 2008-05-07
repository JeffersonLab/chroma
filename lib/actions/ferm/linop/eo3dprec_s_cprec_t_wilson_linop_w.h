#ifndef EO3DPREC_S_CPREC_T_WILSON_LINOP_W_H
#define EO3DPREC_S_CPREC_T_WILSON_LINOP_W_H

#include "qdp_config.h"

#if QDP_ND == 4
#if QDP_NC == 3
#if QDP_NS == 4

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
  
  class EO3DPrecSCprecTWilsonLinOp : 
    public EO3DPrecSpaceCentralPrecTimeLinearOperator<LatticeFermion, 
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
    EO3DPrecSCprecTWilsonLinOp() {}


    //! Full constructor with Anisotropy
    EO3DPrecSCprecTWilsonLinOp(Handle< FermState<T,P,Q> > fs_,
		      const Real& Mass_,
		      const AnisoParam_t& aniso_)
      {create(fs_,Mass_,aniso_);}

    //! Destructor is automatic
    ~EO3DPrecSCprecTWilsonLinOp() {}

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
    void invCLeftLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const;

    //! Apply inv (C_R)^{-1}
    void invCRightLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const;

    //! Apply C_L
    void cLeftLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const;

    //! Apply C_R
    void cRightLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const;

    inline
    void evenEvenLinOp(T& chi, const T& psi, enum PlusMinus isign) const { 
      chi[rb3[0]] = psi;
    }

    inline
    void evenEvenInvLinOp(T& chi, const T& psi, enum PlusMinus isign) const { 
      chi[rb3[0]] = psi;
    }

    inline
    void oddOddLinOp(T& chi, const T& psi, enum PlusMinus isign) const { 
      chi[rb3[1]] = psi;
    }

    inline 
    void evenOddLinOp(T& chi, const T& psi, enum PlusMinus isign) const { 
      T tmp1, tmp2;
      
      switch(isign) { 
      case PLUS:
	{
	
	  cRightLinOp(tmp1, psi, isign, 1);
	  spaceLinOp(tmp2, tmp1, isign, 0);
	  cLeftLinOp(chi, tmp2, isign, 0);
	}
	break;
      case MINUS:
	{
	  cLeftLinOp(tmp1, psi, isign, 1);
	  spaceLinOp(tmp2, tmp1, isign, 0);
	  cRightLinOp(chi, tmp2, isign, 0);
	}
	break;
      };
    }


    inline 
    void oddEvenLinOp(T& chi, const T& psi, enum PlusMinus isign) const { 
      T tmp1, tmp2;
      
      switch(isign) { 
      case PLUS:
	{
	
	  cRightLinOp(tmp1, psi, isign, 0);
	  spaceLinOp(tmp2, tmp1, isign, 1);
	  cLeftLinOp(chi, tmp2, isign, 1);
	}
	break;
      case MINUS:
	{
	  cLeftLinOp(tmp1, psi, isign, 0);
	  spaceLinOp(tmp2, tmp1, isign, 1);
	  cRightLinOp(chi, tmp2, isign, 1);
	}
	break;
      };
    }


    // Override: Get rid of the inv calls
    
    //! Apply the operator onto a source vector
    virtual void operator() (T& chi, const T& psi, enum PlusMinus isign) const
    {
      T   tmp1, tmp2; moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      switch (isign)
      {
      case PLUS:
	// \chi = ( M_oo - M_oe M^{-1}_ee M_eo ) \psi
	
	// chi = M_oo \psi
	oddOddLinOp(chi, psi, PLUS);

	// tmp2 = M_oe M^{-1}_ee M_eo  \psi
	evenOddLinOp(tmp1, psi, PLUS);
	oddEvenLinOp(tmp2, tmp1, PLUS);
	
	chi[rb3[1]] -= tmp2;
	break;

      case MINUS:
	// \chi = ( M_oo - M_oe M^{-1}_ee M_eo )^\dagger \psi
	//      = M^\dagger_oo \psi - M^\dagger_{oe} ( M^{-1}_ee )^\dagger M^\dagger{eo}	
	//
	// NB: Daggering acts on checkerboarding to achieve the result above.

	oddOddLinOp(chi, psi, MINUS);

	// tmp2 = M_oe M^{-1}_ee M_eo  \psi
	evenOddLinOp(tmp1, psi, MINUS);
	oddEvenLinOp(tmp2, tmp1, MINUS);

	chi[rb3[1]] -= tmp2;
	break;

      default:
	QDPIO::cerr << "unknown sign" << endl;
	QDP_abort(1);
      }

      getFermBC().modifyF(chi, QDP::rb3[1]);
    }


    //! Apply the d/dt of the preconditioned linop
    void deriv(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const {

      QDPIO::cerr << "Not Yet Implemented" << endl;
      QDP_abort(1);

    }
    
    //! Get log det ( T^\dag T )
    Double logDetTDagT(void) const {
      QDPIO::cerr << "Not Yet Implemented" << endl;
      QDP_abort(1);
      return (double)0;
    }


    //! Get the force due to the det T^\dag T bit
    void derivLogDetTDagT(P& ds_u, enum PlusMinus isign) const {
      QDPIO::cerr << "Not Yet Implemented" << endl;
      QDP_abort(1);
    }


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
      getFermBC().modifyF(chi, rb3[cb3d]);
    }

    int getTMax() const { 
      return Layout::lattSize()[Nd-1];
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
  };

} // End Namespace Chroma


#endif
#endif
#endif

#endif
