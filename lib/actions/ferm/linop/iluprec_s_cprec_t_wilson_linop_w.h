#ifndef ILUPREC_S_CPREC_T_WILSON_LINOP_W_H
#define ILUPREC_S_CPREC_T_WILSON_LINOP_W_H

#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3

#include "linearop.h"
#include "central_tprec_linop.h"

#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/iluprec_s_cprec_t_wilsonlike_linop_w.h"

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
    public ILUPrecSCprecTWilsonLikeLinOp 
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


    //! Apply the the space block onto a source vector
    //  cb3d is the 3d (rb3) checkerboard of the target
    void spaceLinOp(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const {
      Real mhalf=Real(-0.5);
      Dw3D.apply(chi, psi, isign, cb3d);
      chi[ rb3[cb3d] ] *= mhalf;
      getFermBC().modifyF(chi);

    }

    void derivSpaceOp(P& ds_u, const T& X, const T& Y, enum PlusMinus isign, int cb3d) const 
    {
      Real mhalf=Real(-0.5);
      Dw3D.deriv(ds_u, X, Y, isign, cb3d);
      for(int mu=0; mu < 3; mu++) { 
	ds_u[mu]*= mhalf;
      }
      ds_u[3]=zero;
      getFermBC().zero(ds_u);

    }

   // A = 0 so A-1 = -1 => (A-1) psi = -psi and also for the dagger.
    inline
    void AMinusOneOp(T& chi, const T& psi, enum PlusMinus isign) const { 
      chi = Real(-1)*psi;
      getFermBC().modifyF(chi);
    }
 
    inline
    void derivAMinusOne(P& ds_u, const T& X, const T& Y, enum PlusMinus isign) const 
    {
      // Trivial for Wilson Case
      ds_u.resize(Nd);
      for(int mu = 0; mu < Nd; mu++) { 
	ds_u[mu] = zero;
      }
      getFermBC().zero(ds_u);
    }
    
    //! Get log det ( T^\dag T )
    Double logDetTDagT(void) const {
      return logDetTSq;
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
