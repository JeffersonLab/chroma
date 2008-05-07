#ifndef EO3DPREC_S_CPREC_T_CLOVER_LINOP_W_H
#define EO3DPREC_S_CPREC_T_CLOVER_LINOP_W_H

#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3

#include "linearop.h"
#include "central_tprec_linop.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/central_tprec_nospin_utils.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "actions/ferm/invert/invcg2.h"
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
  class EELinOp : public LinearOperator<LatticeFermion> {
  public:
    ~EELinOp() {}
    EELinOp(const EO3DPrecSpaceCentralPrecTimeLinearOperator<LatticeFermion, 
	      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& EO_,
	    enum PlusMinus isign) : EOLinOp(EO_), my_isign(isign) {}
    
    virtual void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const {
      switch( my_isign ) { 
      case PLUS:
	EOLinOp.evenEvenLinOp(chi, psi, isign);
	break;
      case MINUS:
	{
	  enum PlusMinus opp_sign = (isign == PLUS ? MINUS : PLUS);
	  EOLinOp.evenEvenLinOp(chi, psi, opp_sign);
	}
	break;
      default:
	break;
      };
    }
    //! Return the subset on which the operator acts
    virtual const Subset& subset() const { 
      return rb3[0];
    }


  private:
    const EO3DPrecSpaceCentralPrecTimeLinearOperator<LatticeFermion, 
      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& EOLinOp;

    enum PlusMinus my_isign;
  };
 
  class EO3DPrecSCprecTCloverLinOp : 
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
    EO3DPrecSCprecTCloverLinOp() {}


    //! Full constructor with Anisotropy
    EO3DPrecSCprecTCloverLinOp(Handle< FermState<T,P,Q> > fs_,
			       const CloverFermActParams& param_)
      {create(fs_,param_);}

    //! Destructor is automatic
    ~EO3DPrecSCprecTCloverLinOp() {}

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

    int getTMax() const { 
      return Layout::lattSize()[tDir()];
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
    void diagLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3) const { 
      // This is now  1 + C_L (A-fact) C_R on even even
      T tmp1 = zero;
      T tmp2 = zero; 

      switch(isign) { 
      case PLUS: 
	{
	  cRightLinOp(tmp1, psi, PLUS, cb3);

	  const int* tab = rb3[cb3].siteTable().slice();
	  for(int j=0; j < rb3[cb3].siteTable().size(); j++) { 
	    int site=tab[j];
	    APlusFact.applySite(tmp2, tmp1, PLUS, site);
	  }
	  tmp2[rb3[cb3]] -= fact * tmp1;

	  cLeftLinOp(tmp1, tmp2, PLUS, cb3);

	  chi[rb3[cb3]] = psi + tmp1;
	}
	break;
      case MINUS:
	{
	  cLeftLinOp(tmp1, psi, MINUS, cb3);
	  const int* tab = rb3[cb3].siteTable().slice();
	  for(int j=0; j < rb3[cb3].siteTable().size(); j++) { 
	    int site=tab[j];
	    APlusFact.applySite(tmp2, tmp1, MINUS, site);
	  }

	  tmp2[rb3[cb3]] -= fact * tmp1;
	  cRightLinOp(tmp1, tmp2, MINUS, cb3);
	  chi[rb3[cb3]] = psi + tmp1;
	}
	break;
      default:
	QDPIO::cout << "error: Unknown ISIGN" << endl;
	QDP_abort(1);
	break;
      }
    }



    inline
    void evenEvenLinOp(T& chi, const T& psi, enum PlusMinus isign) const { 
      diagLinOp(chi,psi, isign, 0);
    }

    inline
    void evenEvenInvLinOp(T& chi, const T& psi, enum PlusMinus isign) const { 
      // Do with a CG?
      EELinOp eeLin(*this, isign);
      enum PlusMinus opp_sign = (isign == PLUS ? MINUS : PLUS );

      T tmp1;
      chi[rb3[0]]=zero;
      evenEvenLinOp(tmp1, psi, opp_sign);
      Real RsdCG=Real(1.0e-8);
      int MaxCG=200;
      InvCG2( eeLin, tmp1, chi, RsdCG, MaxCG);

    }

    inline
    void oddOddLinOp(T& chi, const T& psi, enum PlusMinus isign) const { 
      diagLinOp(chi,psi, isign, 1);
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

    unsigned long nFlops() const 
    { 
      //      QDPIO::cout << "Flopcount Not Yet Implemented " << endl;
      return 0;
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

    
  private:
    
    //! Apply the the space block onto a source vector
    //  cb3d is the 3d (rb3) checkerboard of the target
    void spaceLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const {
      Real mhalf=Real(-0.5);
      Dw3D.apply(chi, psi, isign, cb3d);
      chi[rb3[cb3d]] *= mhalf;
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

    CloverTerm APlusFact;
    CloverFermActParams param;

    WilsonDslash3D Dw3D;

  };

} // End Namespace Chroma

#endif
#endif
#endif

#endif
