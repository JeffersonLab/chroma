#ifndef ILU2PREC_S_CPREC_T_WILSON_LINOP_W_H
#define ILU2PREC_S_CPREC_T_WILSON_LINOP_W_H

#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3

#include "linearop.h"
#include "central_tprec_linop.h"

#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/ilu2prec_s_cprec_t_wilsonlike_linop_w.h"

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
  
  class ILU2PrecSCprecTWilsonLinOp : 
    public ILU2PrecSCprecTWilsonLikeLinOp 
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;
    typedef PScalar<PColorMatrix<RComplex<REAL>, 3> > CMat;              // Useful type: ColorMat with no Outer<>
    typedef PSpinVector<PColorVector<RComplex<REAL>, 3>, 2> HVec_site;   // Useful type: Half Vec with no Outer<>
    //! Partial constructor
    ILU2PrecSCprecTWilsonLinOp() {}


    //! Full constructor with Anisotropy
    ILU2PrecSCprecTWilsonLinOp(Handle< FermState<T,P,Q> > fs_,
		      const Real& Mass_,
		      const AnisoParam_t& aniso_)
      {create(fs_,Mass_,aniso_);}

    //! Destructor is automatic
    ~ILU2PrecSCprecTWilsonLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return fs->getBC();}


    //! Creation routine with Anisotropy
    void create(Handle< FermState<T,P,Q> > fs_,
		const Real& Mass_,
		const AnisoParam_t& aniso_);

#if 0
    void operator() (T& chi, const T& psi, enum PlusMinus isign) const
    {
      QDPIO::cout << "Foo" << endl;
      T   tmp1, tmp2,tmp3; moveToFastMemoryHint(tmp1); moveToFastMemoryHint(tmp2);

      switch (isign)
      {
      case PLUS:
	chi = psi;
	DBar(tmp1, psi,  PLUS, 0);
	DBar(tmp2, tmp1, PLUS, 1);
	chi[rb3[1]] -= tmp2;

	break;

	
      case MINUS:
	chi = psi;
	DBar(tmp1, psi,  MINUS, 0);
	DBar(tmp2, tmp1, MINUS, 1);
	chi[rb3[1]] -= tmp2;
	
	break;

      default:
	QDPIO::cerr << "unknown sign" << endl;
	QDP_abort(1);
      }
      getFermBC().modifyF(chi);

    }
#endif
    //! Apply the the space block onto a source vector
    //  cb3d is the 3d (rb3) checkerboard of the target
    void Dslash3D(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const {
      Real mhalf=Real(0.5); // Minus sign explicit in the preconditioner. But factor of 2 
                            // Still here...
      Dw3D.apply(chi, psi, isign, cb3d);
      chi[ rb3[cb3d] ] *= mhalf;
      getFermBC().modifyF(chi);

    }

    inline
      void ABar(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const { 
      chi[rb3[cb3d]] = zero;
      getFermBC().modifyF(chi);
    }

    inline
      void AH(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const { 
      chi[rb3[cb3d]] = zero;
      getFermBC().modifyF(chi);
    }

 
    const multi3d< int >&  getTsiteArray(void) const 
    {
      return tsite;
    }

    const multi3d< CMat >& getPMatrixArray(void) const
    {
      return P_mat;
    }

    const multi3d< CMat >& getPMatrixDaggerArray(void) const 
    {
      return P_mat_dag;
    }

    const multi2d< CMat >& getQMatrixInvArray(void) const 
    {
      return Q_mat_inv;
    }

    const multi2d< CMat >& getQMatrixDaggerInvArray(void) const 
    {
      return Q_mat_dag_inv;
    }

    const Real& getFactor() const 
    {
      return fact;
    }

    const Real& getInvFactor() const 
    {
      return invfact;
    }

    const multi1d<LatticeColorMatrix>& getLinks(void) const { 
      return u;
    }

    inline bool schroedingerTP() const { 
      return schrTP;
    }

  protected:
    inline int getTMax(void) const {
      // Always return Nt for now
      return Layout::lattSize()[Nd-1];
    }
  private:
    bool schrTP;
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
