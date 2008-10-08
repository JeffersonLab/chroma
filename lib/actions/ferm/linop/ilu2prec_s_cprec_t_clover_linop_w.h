#ifndef ILU2PREC_S_CPREC_T_CLOVER_LINOP_W_H
#define ILU2PREC_S_CPREC_T_CLOVER_LINOP_W_H
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
#include "actions/ferm/linop/ilu2prec_s_cprec_t_wilsonlike_linop_w.h"

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
  
  class ILU2PrecSCprecTCloverLinOp : public ILU2PrecSCprecTWilsonLikeLinOp
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;
    typedef PScalar<PColorMatrix<RComplex<REAL>, 3> > CMat;              // Useful type: ColorMat with no Outer<>
    typedef PSpinVector<PColorVector<RComplex<REAL>, 3>, 2> HVec_site;   // Useful type: Half Vec with no Outer<>
    //! Partial constructor
    ILU2PrecSCprecTCloverLinOp() {}


    //! Full constructor with Anisotropy
    ILU2PrecSCprecTCloverLinOp(Handle< FermState<T,P,Q> > fs_,
			      const CloverFermActParams& param_ )
      {create(fs_,param_);}

    //! Destructor is automatic
    ~ILU2PrecSCprecTCloverLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return fs->getBC();}


    //! Creation routine with Anisotropy
    void create(Handle< FermState<T,P,Q> > fs_,
		const CloverFermActParams& param_);


  protected:
    
    //! Apply the the space block onto a source vector
    //  cb3d is the 3d (rb3) checkerboard of the target
    void Dslash3D(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const {
      Real mhalf=Real(0.5); // - sign is now in preconditioners
      Dw3D.apply(chi, psi, isign, cb3d);
      chi[ rb3[cb3d] ] *= mhalf;
      getFermBC().modifyF(chi, rb3[cb3d]);
    }

    //! Do chi = (A - 1) psi with A = -(c_sw/4) sigma_munu F_munu
    inline
      void AH(T& chi, const T& psi, enum PlusMinus isign, int cb3d) const { 
      T tmp1;

      // The clover term does A_plus_fact = factI - (c_sw/4) sigma_munu F_munu = A + fact
      // where the desired A is just A = - (c_sw/4) sigma_munu F_munu

      // I need to get rid of the factI. Luckily, since I already
      // have to do A - 1, I can just change this to A-(1 + fact)I

      // tmp = (A + fact) psi
      // Yucky...
      const int *tab = rb3[cb3d].siteTable().slice();
      for(int j=0; j < rb3[cb3d].siteTable().size(); j++) { 
	int site = tab[j];
	APlusFact.applySite(tmp1, psi, isign, site);
      }

      // chi
      chi[rb3[cb3d]] = tmp1 - fact*psi;
      getFermBC().modifyF(chi);
    }
    
    //! Flopcounter
    unsigned long nFlops() const 
    { 
      //      QDPIO::cout << "Flopcount Not Yet Implemented " << endl;
      return 0;
    }
    
    // Fill out the blanks
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

    inline int getTMax(void) const {
      return t_max;
    }
  private:
    bool schrTP;

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
    Double logDetTSq;
    int t_max;
  };

} // End Namespace Chroma

#endif
#endif
#endif

#endif
