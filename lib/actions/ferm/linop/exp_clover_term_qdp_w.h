// -*- C++ -*-
/*! \file
 *  \brief Clover term linear operator
 */

#ifndef __exp_clover_term_qdp_w_h__
#define __exp_clover_term_qdp_w_h__

#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/exp_clover_term_base_w.h"
#include "meas/glue/mesfield.h"
#include "qdp_allocator.h"
#include "state.h"
#include <cmath>
#include <complex>
namespace Chroma
{

  namespace
  {

    template <typename T>
    struct ExpClovTriang {

      // Unexponentiated part
      PrimitiveClovTriang<T> A;

      /*! Linear combination coefficients to generate the exponential */
      RScalar<T> q[2][6];

      // Exponentiated part
      RScalar<T> qinv[2][6];

      /*! Coefficients for the force*/
      RScalar<T> C[2][6][6];
    };

    /*! This accessor class allows me a convenient way to acces the
        diagonal + lower diagonal storage for the hermitian matrix */
    template <typename T, int block>
    struct ClovAccessor {
      ClovAccessor(PrimitiveClovTriang<T>& t) : tri(t)
      {
      }

      inline RComplex<T> operator()(int row, int col) const
      {
	RComplex<T> ret_val;
	if (row == col)
	{
	  // Diagonal Piece:
	  ret_val = tri.diag[block][row];
	}
	else if (row > col)
	{
	  // Lower triangular portion
	  ret_val = tri.offd[block][(row * (row - 1)) / 2 + col];
	}
	else if (row < col)
	{
	  // Upper triangular portion: transpose ( row <-> col) and conjugate
	  ret_val = conj(tri.offd[block][(col * (col - 1)) / 2 + row]);
	}

	return ret_val;
      }

      inline void insert(int row, int col, const RComplex<T>& value)
      {
	if (row == col)
	{
	  // Diagonal piece -- must be real.
	  tri.diag[block][row] = RScalar<T>(real(value));
	}
	else if (row > col)
	{
	  // Lower triangular portion
	  tri.offd[block][(row * (row - 1)) / 2 + col] = value;
	}
	else if (row < col)
	{
	  // Upper triangular portion: transpose ( row <-> col) and conjugate
	  tri.offd[block][(col * (col - 1)) / 2 + row] = conj(value);
	}
      }

    private:
      // A reference to the triangular storage
      PrimitiveClovTriang<T>& tri;
    };

    template <typename T, int block>
    struct Traces {
      Traces(ExpClovTriang<T>& E_) : E(E_)
      {
      }

      // Simple mat mult routine
      inline void multiply(ClovAccessor<T, block>& out, const ClovAccessor<T, block>& M1,
			   const ClovAccessor<T, block>& M2)
      {
	RComplex<T> zip(RScalar<T>((T)0), RScalar<T>((T)0));

	// NB: We only need to compute the diagonal and lower diagonal
	// elements because the matrices are hermitiean.
	for (int row = 0; row < 6; ++row)
	{
	  for (int col = 0; col <= row; ++col)
	  {
	    // Pour row down column
	    RComplex<T> dotprod = zip;
	    for (int k = 0; k < 6; ++k)
	    {
	      dotprod += M1(row, k) * M2(k, col);
	    }
	    out.insert(row, col, dotprod);
	  }
	}
      }

      // Simple mat mult routine
      inline void copy(ClovAccessor<T, block>& out, const ClovAccessor<T, block>& in)
      {
	// NB: We only need to compute the diagonal and lower diagonal
	// elements because the matrices are hermitiean.
	for (int row = 0; row < 6; ++row)
	{
	  for (int col = 0; col <= row; ++col)
	  {
	    out.insert(row, col, in(row, col));
	  }
	}
      }

      inline void traces(multi1d<RScalar<T>>& tr)
      {
	RComplex<T> zip(RScalar<T>((T)0), RScalar<T>((T)0));

	// The first 5 will map onto the hither powers.
	ClovAccessor<T, block> A(E.A);

	// Trace A
	tr.resize(6);
	tr[0] = RScalar<T>((T)0);
	for (int i = 0; i < 6; i++)
	{
	  tr[0] += real(A(i, i));
	}

	PrimitiveClovTriang<T> prev;
	PrimitiveClovTriang<T> curr;
	ClovAccessor<T, block> Prev(prev);
	ClovAccessor<T, block> Curr(curr);

	copy(Prev, A);

	for (int pow = 1; pow <= 5; pow++)
	{
	  multiply(Curr, Prev, A);
	  tr[pow] = RScalar<T>((T)0);
	  for (int i = 0; i < 6; i++)
	  {
	    tr[pow] += real(Curr(i, i));
	  }
	  copy(Prev, Curr);
	}
      }

    private:
      ExpClovTriang<T>& E;
    };

    constexpr int N_exp_default = 17;

    template <typename REALT, int block>
    inline void siteApplicationBlock(RComplex<REALT>* __restrict__ cchi,
				     const PrimitiveClovTriang<REALT>& tri_in,
				     const RComplex<REALT>* const __restrict__ ppsi)
    {

      if constexpr (block == 0)
      {
	cchi[0] = tri_in.diag[0][0] * ppsi[0] + conj(tri_in.offd[0][0]) * ppsi[1] +
		  conj(tri_in.offd[0][1]) * ppsi[2] + conj(tri_in.offd[0][3]) * ppsi[3] +
		  conj(tri_in.offd[0][6]) * ppsi[4] + conj(tri_in.offd[0][10]) * ppsi[5];

	cchi[1] = tri_in.diag[0][1] * ppsi[1] + tri_in.offd[0][0] * ppsi[0] +
		  conj(tri_in.offd[0][2]) * ppsi[2] + conj(tri_in.offd[0][4]) * ppsi[3] +
		  conj(tri_in.offd[0][7]) * ppsi[4] + conj(tri_in.offd[0][11]) * ppsi[5];

	cchi[2] = tri_in.diag[0][2] * ppsi[2] + tri_in.offd[0][1] * ppsi[0] +
		  tri_in.offd[0][2] * ppsi[1] + conj(tri_in.offd[0][5]) * ppsi[3] +
		  conj(tri_in.offd[0][8]) * ppsi[4] + conj(tri_in.offd[0][12]) * ppsi[5];

	cchi[3] = tri_in.diag[0][3] * ppsi[3] + tri_in.offd[0][3] * ppsi[0] +
		  tri_in.offd[0][4] * ppsi[1] + tri_in.offd[0][5] * ppsi[2] +
		  conj(tri_in.offd[0][9]) * ppsi[4] + conj(tri_in.offd[0][13]) * ppsi[5];

	cchi[4] = tri_in.diag[0][4] * ppsi[4] + tri_in.offd[0][6] * ppsi[0] +
		  tri_in.offd[0][7] * ppsi[1] + tri_in.offd[0][8] * ppsi[2] +
		  tri_in.offd[0][9] * ppsi[3] + conj(tri_in.offd[0][14]) * ppsi[5];

	cchi[5] = tri_in.diag[0][5] * ppsi[5] + tri_in.offd[0][10] * ppsi[0] +
		  tri_in.offd[0][11] * ppsi[1] + tri_in.offd[0][12] * ppsi[2] +
		  tri_in.offd[0][13] * ppsi[3] + tri_in.offd[0][14] * ppsi[4];
      }
      else
      {
	cchi[6] = tri_in.diag[1][0] * ppsi[6] + conj(tri_in.offd[1][0]) * ppsi[7] +
		  conj(tri_in.offd[1][1]) * ppsi[8] + conj(tri_in.offd[1][3]) * ppsi[9] +
		  conj(tri_in.offd[1][6]) * ppsi[10] + conj(tri_in.offd[1][10]) * ppsi[11];

	cchi[7] = tri_in.diag[1][1] * ppsi[7] + tri_in.offd[1][0] * ppsi[6] +
		  conj(tri_in.offd[1][2]) * ppsi[8] + conj(tri_in.offd[1][4]) * ppsi[9] +
		  conj(tri_in.offd[1][7]) * ppsi[10] + conj(tri_in.offd[1][11]) * ppsi[11];

	cchi[8] = tri_in.diag[1][2] * ppsi[8] + tri_in.offd[1][1] * ppsi[6] +
		  tri_in.offd[1][2] * ppsi[7] + conj(tri_in.offd[1][5]) * ppsi[9] +
		  conj(tri_in.offd[1][8]) * ppsi[10] + conj(tri_in.offd[1][12]) * ppsi[11];

	cchi[9] = tri_in.diag[1][3] * ppsi[9] + tri_in.offd[1][3] * ppsi[6] +
		  tri_in.offd[1][4] * ppsi[7] + tri_in.offd[1][5] * ppsi[8] +
		  conj(tri_in.offd[1][9]) * ppsi[10] + conj(tri_in.offd[1][13]) * ppsi[11];

	cchi[10] = tri_in.diag[1][4] * ppsi[10] + tri_in.offd[1][6] * ppsi[6] +
		   tri_in.offd[1][7] * ppsi[7] + tri_in.offd[1][8] * ppsi[8] +
		   tri_in.offd[1][9] * ppsi[9] + conj(tri_in.offd[1][14]) * ppsi[11];

	cchi[11] = tri_in.diag[1][5] * ppsi[11] + tri_in.offd[1][10] * ppsi[6] +
		   tri_in.offd[1][11] * ppsi[7] + tri_in.offd[1][12] * ppsi[8] +
		   tri_in.offd[1][13] * ppsi[9] + tri_in.offd[1][14] * ppsi[10];
      }
    }

#if 0

    template<typename T>
    inline
    void siteExponentiate(ExpClovTriang<T>& tri_in)
    {
      RComplex<T> zip( RScalar<T>((T)0), RScalar<T>((T)0) );
      for(int block=0; block < 2; ++block) {
        // q0 * I -- no offdiag piece in I

        for(int i=0; i < 6; ++i) {
          tri_in.Exp[0].diag[block][i] = tri_in.q[block][0];
          tri_in.Exp[1].diag[block][i] = tri_in.qinv[block][0];
        }

        for(int ord=1; ord <= 5; ++ord) {
          for(int i=0; i < 6; ++i) {
            tri_in.Exp[0].diag[block][i] += tri_in.q[block][ord]*tri_in.A[ord-1].diag[block][i];
            tri_in.Exp[1].diag[block][i] += tri_in.qinv[block][ord]*tri_in.A[ord-1].diag[block][i];
          }
        }

        for(int i=0; i < 15; ++i) {
          tri_in.Exp[0].offd[block][i] = zip;
          tri_in.Exp[1].offd[block][i] = zip;
        }

        for(int ord=1; ord <= 5; ++ord) {
          for(int i=0; i < 15; ++i) {
            tri_in.Exp[0].offd[block][i] += tri_in.q[block][ord]*tri_in.A[ord-1].offd[block][i];
            tri_in.Exp[1].offd[block][i] += tri_in.qinv[block][ord]*tri_in.A[ord-1].offd[block][i];
          }
        }

       
      }
    }
#endif

    template <typename REALT, int i = 0>
    inline void siteApplicationExp(RComplex<REALT>* __restrict__ cchi,
				   const ExpClovTriang<REALT>& tri_in,
				   const RComplex<REALT>* const __restrict__ ppsi)
    {

#if 0
      siteApplicationBlock<REALT,0>(cchi,tri_in.Exp[i],ppsi);
      siteApplicationBlock<REALT,1>(cchi,tri_in.Exp[i],ppsi);
#endif

#if 0
     // Accumulate exponential from stored powers of A
     RComplex<REALT> tmp[12];
      for(int i=0; i < 6; ++i) {
        cchi[i] = tri_in.q[0][0]*ppsi[i];
      }
      for(int pow=1; pow <=5; pow++) {
        siteApplicationBlock<REAL,0>(tmp,tri_in.A[pow-1],ppsi);
        for(int i=0;i < 6; ++i) {
          cchi[i] += tri_in.q[0][pow]*tmp[i];
        }
      }

        
      for(int i=6; i < 12; ++i) {
        cchi[i] = tri_in.q[1][0]*ppsi[i];
      }
      for(int pow=1; pow <=5; pow++) {
        siteApplicationBlock<REAL,1>(tmp,tri_in.A[pow-1],ppsi);
        for(int i=6; i < 12; ++i) {
          cchi[i] += tri_in.q[1][pow]*tmp[i];
        }
      }
#endif

#if 1

      // Accumulate exponential from only A
      RComplex<REALT> tmp[12];
      // Top block
      // chi = psi
      for (int cspin = 0; cspin < 6; ++cspin)
      {
	cchi[cspin] = ppsi[cspin];
      }

      // Main loop:  chi = psi + q[i]/q[i-1] A chi
      for (int pow = 5; pow > 0; --pow)
      {
	siteApplicationBlock<REALT, 0>(tmp, tri_in.A, cchi);
	for (int cspin = 0; cspin < 6; cspin++)
	{

	  if constexpr (i == 0)
	  {
	    // Operator
	    cchi[cspin] = ppsi[cspin] + (tri_in.q[0][pow] / tri_in.q[0][pow - 1]) * tmp[cspin];
	  }
	  else
	  {
	    // Inverse
	    cchi[cspin] =
	      ppsi[cspin] + (tri_in.qinv[0][pow] / tri_in.qinv[0][pow - 1]) * tmp[cspin];
	  }
	}
      }

      for (int cspin = 0; cspin < 6; ++cspin)
      {
	if constexpr (i == 0)
	{
	  // Operator
	  cchi[cspin] *= tri_in.q[0][0];
	}
	else
	{
	  // inverse
	  cchi[cspin] *= tri_in.qinv[0][0];
	}
      }

      // Second Block
      // chi = psi
      for (int cspin = 6; cspin < 12; ++cspin)
      {
	cchi[cspin] = ppsi[cspin];
      }

      // Main loop:  chi = psi + q[i]/q[i-1] A chi
      for (int pow = 5; pow > 0; --pow)
      {
	siteApplicationBlock<REALT, 1>(tmp, tri_in.A, cchi);

	for (int cspin = 6; cspin < 12; cspin++)
	{
	  if constexpr (i == 0)
	  {
	    // Operator
	    cchi[cspin] = ppsi[cspin] + (tri_in.q[1][pow] / tri_in.q[1][pow - 1]) * tmp[cspin];
	  }
	  else
	  {
	    // Inverse
	    cchi[cspin] =
	      ppsi[cspin] + (tri_in.qinv[1][pow] / tri_in.qinv[1][pow - 1]) * tmp[cspin];
	  }
	}
      }

      for (int cspin = 6; cspin < 12; ++cspin)
      {

	if constexpr (i == 0)
	{
	  // Operator
	  cchi[cspin] *= tri_in.q[1][0];
	}
	else
	{
	  // Inverse
	  cchi[cspin] *= tri_in.qinv[1][0];
	}
      }
#endif
    }

    template <typename REALT>
    inline void siteApplicationPower(RComplex<REALT>* __restrict__ cchi,
				     const ExpClovTriang<REALT>& tri_in,
				     const RComplex<REALT>* const __restrict__ ppsi, int power)
    {
      if (power == 0)
      {
	// Straight Copy
	for (int i = 0; i < 12; ++i)
	  cchi[i] = ppsi[i];
      }
      else
      {
	// Higher powers, just construct as before.
	RComplex<REALT> tmp[12];

	siteApplicationBlock<REALT, 0>(cchi, tri_in.A, ppsi);
	for (int p = power; p > 1; --p)
	{
	  for (int i = 0; i < 6; ++i)
	    tmp[i] = cchi[i];
	  siteApplicationBlock<REALT, 0>(cchi, tri_in.A, tmp);
	}

	siteApplicationBlock<REALT, 1>(cchi, tri_in.A, ppsi);
	for (int p = power; p > 1; --p)
	{
	  for (int i = 6; i < 12; ++i)
	    tmp[i] = cchi[i];
	  siteApplicationBlock<REALT, 1>(cchi, tri_in.A, tmp);
	}
      }
    }

  }
  // Reader/writers
  /*! \ingroup linop */

  //! Clover term
  /*!
   * \ingroup linop
   *
   */
  template <typename T, typename U, int N_exp = N_exp_default>
  class QDPExpCloverTermT : public ExpCloverTermBase<T, U>
  {
  public:
    // Typedefs to save typing
    typedef typename WordType<T>::Type_t REALT;

    typedef OLattice<PScalar<PScalar<RScalar<typename WordType<T>::Type_t>>>> LatticeREAL;
    typedef OScalar<PScalar<PScalar<RScalar<REALT>>>> RealT;

    //! Empty constructor. Must use create later
    QDPExpCloverTermT();

    //! No real need for cleanup here
    ~QDPExpCloverTermT()
    {
      if (tri != nullptr)
      {
	QDP::Allocator::theQDPAllocator::Instance().free(tri);
      }
    }

    //! Creation routine
    void create(Handle<FermState<T, multi1d<U>, multi1d<U>>> fs, const CloverFermActParams& param_);

    virtual void create(Handle<FermState<T, multi1d<U>, multi1d<U>>> fs,
			const CloverFermActParams& param_, const QDPExpCloverTermT<T, U>& from_);

    //! Computes the inverse of the term on cb using Cholesky
    /*!
     * \param cb   checkerboard of work (Read)
     */
    void choles(int cb) override;

    //! Computes the inverse of the term on cb using Cholesky
    /*!
     * \param cb   checkerboard of work (Read)
     * \return logarithm of the determinant  
     */
    Double cholesDet(int cb) const override;
    /**
     * Apply a dslash
     *
     * Performs the operation
     *
     *  chi <-   (L + D + L^dag) . psi
     *
     * where
     *   L       is a lower triangular matrix
     *   D       is the real diagonal. (stored together in type TRIANG)
     *
     * Arguments:
     * \param chi     result                                      (Write)
     * \param psi     source                                      (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
     */

    void fillRefDiag(Real diag);

    // Reference exponential using old fashioned taylor expansion
    void applyRef(T& chi, const T& psi, enum PlusMinus isign, int N = N_exp) const;

    // Appl;y a power of a matrix from A^0 to A^5
    void applyPowerSite(T& chi, const T& psi, enum PlusMinus isign, int site, int power = 1) const;

    // Appl;y a power of a matrix from A^0 to A^5
    void applyPower(T& chi, const T& psi, enum PlusMinus isign, int cb, int power = 1) const;

    // Apply exponential operator
    void apply(T& chi, const T& psi, enum PlusMinus isign, int cb) const override;

    // Apply exponential operator
    void applyInv(T& chi, const T& psi, enum PlusMinus isign, int cb) const;

    // Apply exponential operator at a site
    void applySite(T& chi, const T& psi, enum PlusMinus isign, int site) const override;

    // Apply exponential operator at a site
    void applySiteInv(T& chi, const T& psi, enum PlusMinus isign, int site) const;

    inline void applyUnexp(T& chi, const T& psi, enum PlusMinus isign, int cb) const
    {
      applyPower(chi, psi, isign, cb, 1); // Explicily apply just the clover term.
    }

    //! Calculates Tr_D ( Gamma_mat L )
    void triacntr(U& B, int mat, int cb) const override;

    //! Return the fermion BC object for this linear operator
    const FermBC<T, multi1d<U>, multi1d<U>>& getFermBC() const override
    {
      return *fbc;
    }

    //! PACK UP the Clover term for QUDA library:
    void packForQUDA(multi1d<QUDAPackedClovSite<REALT>>& quda_pack, int cb) const;
    void tracePowers();

  protected:
    //! Create the clover term on cb
    /*!
     *  \param f         field strength tensor F(mu,nu)        (Read)
     *  \param cb        checkerboard                          (Read)
     */
    void makeClov(const multi1d<U>& f, const RealT& diag_mass);

    //void tracePowers();

    //! Compute the approximation coefficient
    void exponentiate();

    //! Invert the clover term on cb
    void chlclovms(LatticeREAL& tr_Minv, int cb);
    void ldagdlinv(LatticeREAL& tr_Minv, int cb);

    //! Get the u field
    const multi1d<U>& getU() const override
    {
      return u;
    }

    //! Calculates Tr_D ( Gamma_mat L )
    Real getCloverCoeff(int mu, int nu) const override;

  private:
    Handle<FermBC<T, multi1d<U>, multi1d<U>>> fbc;
    multi1d<U> u;
    CloverFermActParams param;
    LatticeDouble tr_M; // Fill this out during create

    ExpClovTriang<REALT>* tri;
  };

  // Empty constructor. Must use create later
  template <typename T, typename U, int N_exp>
  QDPExpCloverTermT<T, U, N_exp>::QDPExpCloverTermT()
  {
    // Always allocate on construction
    int nodeSites = Layout::sitesOnNode();
    tri = (ExpClovTriang<REALT>*)QDP::Allocator::theQDPAllocator::Instance().allocate(
      nodeSites * sizeof(ExpClovTriang<REALT>), QDP::Allocator::DEFAULT);
  }

  // Now copy
  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::create(Handle<FermState<T, multi1d<U>, multi1d<U>>> fs,
					      const CloverFermActParams& param_,
					      const QDPExpCloverTermT<T, U>& from)
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();
    u.resize(Nd);

    u = fs->getLinks();
    fbc = fs->getFermBC();
    param = param_;

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "QDPCloverTerm: error: fbc is null" << std::endl;
      QDP_abort(1);
    }

    //
    // Yuk. Some bits of knowledge of the dslash term are buried in the
    // effective mass term. They show up here. If I wanted some more
    // complicated dslash then this will have to be fixed/adjusted.
    //
    RealT diag_mass;
    {
      RealT ff = param.anisoParam.anisoP ? param.anisoParam.nu / param.anisoParam.xi_0 : Real(1);
      diag_mass = 1 + (Nd - 1) * ff + param.Mass;
    }

    {
      RealT ff = param.anisoParam.anisoP ? Real(1) / param.anisoParam.xi_0 : Real(1);
      param.clovCoeffR *= Real(0.5) * ff / diag_mass;
      param.clovCoeffT *= Real(0.5) / diag_mass;
    }

    /* Calculate F(mu,nu) */
    //multi1d<LatticeColorMatrix> f;
    //mesField(f, u);
    //makeClov(f, diag_mass);

    int nodeSites = Layout::sitesOnNode();
    // Deep copy.
#  pragma omp parallel for
    for (int site = 0; site < nodeSites; ++site)
    {
      tr_M.elem(site).elem().elem().elem() = 0;

      for (int block = 0; block < 2; ++block)
      {

	for (int d = 0; d < 6; ++d)
	{
	  tri[site].A.diag[block][d] = from.tri[site].A.diag[block][d];

	  // Inline accumulate the trace
	  tr_M.elem(site).elem().elem() += tri[site].A.diag[block][d];
	}

	for (int od = 0; od < 15; ++od)
	{
	  tri[site].A.offd[block][od] = from.tri[site].A.offd[block][od];
	}

	// The exponentiation coefficients
	for (int i = 0; i < 6; ++i)
	{
	  tri[site].q[block][i] = from.tri[site].q[block][i];
	}

	for (int i = 0; i < 6; ++i)
	{
	  tri[site].qinv[block][i] = from.tri[site].qinv[block][i];
	}

	// The force coefficients
	for (int i = 0; i < 6; ++i)
	{
	  for (int j = 0; j < 6; ++j)
	  {
	    tri[site].C[block][i][j] = from.tri[site].C[block][i][j];
	  }
	}
      } // End site loop
    }

    END_CODE();
#endif
  }

  //! Creation routine
  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::create(Handle<FermState<T, multi1d<U>, multi1d<U>>> fs,
					      const CloverFermActParams& param_)
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    u.resize(Nd);

    u = fs->getLinks();
    fbc = fs->getFermBC();
    param = param_;

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "QDPCloverTerm: error: fbc is null" << std::endl;
      QDP_abort(1);
    }

    //
    // Yuk. Some bits of knowledge of the dslash term are buried in the
    // effective mass term. They show up here. If I wanted some more
    // complicated dslash then this will have to be fixed/adjusted.
    //
    RealT diag_mass;
    {
      RealT ff = param.anisoParam.anisoP ? param.anisoParam.nu / param.anisoParam.xi_0 : Real(1);
      diag_mass = 1 + (Nd - 1) * ff + param.Mass;
    }

    {
      RealT ff = param.anisoParam.anisoP ? Real(1) / param.anisoParam.xi_0 : Real(1);
      param.clovCoeffR *= RealT(0.5) * ff / diag_mass;
      param.clovCoeffT *= RealT(0.5) / diag_mass;
    }

    /* Calculate F(mu,nu) */
    multi1d<U> f;
    mesField(f, u);
    makeClov(f, diag_mass);
    tracePowers();

#  pragma omp parallel for
    for (int site = 0; site < Layout::sitesOnNode(); ++site)
    {
      tr_M.elem(site).elem().elem().elem() = 0;

      for (int block = 0; block < 2; ++block)
      {
	for (int d = 0; d < 6; ++d)
	{
	  // Inline accumulate the trace
	  tr_M.elem(site).elem().elem() += tri[site].A.diag[block][d];
	}
      }
    }

    END_CODE();
#endif
  }

  namespace QDPExpCloverEnv
  {

    template <typename U>
    struct FillRefArg {
      typedef typename WordType<U>::Type_t REALT;
      typedef OScalar<PScalar<PScalar<RScalar<REALT>>>> RealT;
      const REALT& ref_val;
      ExpClovTriang<REALT>* tri;
    };

    /* This is the extracted site loop for makeClover */
    template <typename U>
    inline void fillRefLoop(int lo, int hi, int myId, FillRefArg<U>* a)
    {
#ifndef QDP_IS_QDPJIT
      typedef typename FillRefArg<U>::RealT RealT;
      typedef typename FillRefArg<U>::REALT REALT;

      const REALT ref_val = a->ref_val;
      ExpClovTriang<REALT>* tri = a->tri;

      // SITE LOOP STARTS HERE
      for (int site = lo; site < hi; ++site)
      {
	for (int chiral = 0; chiral < 2; chiral++)
	{
	  for (int diag_index = 0; diag_index < 2 * Nc; diag_index++)
	  {

	    tri[site].A.diag[chiral][diag_index] = RScalar<REALT>(ref_val);
	  }

	  for (int offdiag_index = 0; offdiag_index < 15; offdiag_index++)
	  {
	    tri[site].A.offd[chiral][offdiag_index].real() = 0;
	    tri[site].A.offd[chiral][offdiag_index].imag() = 0;
	  }
	}
      }
#endif
    }

  } /* end namespace */

  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::fillRefDiag(Real ref_val)
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    if (Ns != 4)
    {
      QDPIO::cerr << __func__ << ": CloverTerm::apply requires Ns==4" << std::endl;
      QDP_abort(1);
    }

    QDPExpCloverEnv::FillRefArg<U> arg = {toDouble(ref_val), tri};
    int num_sites = all.siteTable().size();

    // The dispatch function is at the end of the file
    // ought to work for non-threaded targets too...
    dispatch_to_threads(num_sites, arg, QDPExpCloverEnv::fillRefLoop<U>);
    END_CODE();
#endif
  }

  /*Threaded this. It needs a QMT arg struct and I've extracted the site loop */
  namespace QDPExpCloverEnv
  {

    template <typename U>
    struct QDPCloverMakeClovArg {
      typedef typename WordType<U>::Type_t REALT;
      typedef OScalar<PScalar<PScalar<RScalar<REALT>>>> RealT;
      const RealT& diag_mass;
      const U& f0;
      const U& f1;
      const U& f2;
      const U& f3;
      const U& f4;
      const U& f5;
      ExpClovTriang<REALT>* tri;
    };

    /* This is the extracted site loop for makeClover */
    template <typename U>
    inline void makeClovSiteLoop(int lo, int hi, int myId, QDPCloverMakeClovArg<U>* a)
    {
#ifndef QDP_IS_QDPJIT
      typedef typename QDPCloverMakeClovArg<U>::RealT RealT;
      typedef typename QDPCloverMakeClovArg<U>::REALT REALT;

      const RealT& diag_mass = a->diag_mass;
      const U& f0 = a->f0;
      const U& f1 = a->f1;
      const U& f2 = a->f2;
      const U& f3 = a->f3;
      const U& f4 = a->f4;
      const U& f5 = a->f5;
      ExpClovTriang<REALT>* tri = a->tri;
      const int power = 1;

      // SITE LOOP STARTS HERE
      for (int site = lo; site < hi; ++site)
      {
	for (int jj = 0; jj < 2; jj++)
	{
	  for (int ii = 0; ii < 2 * Nc; ii++)
	  {
	    tri[site].A.diag[jj][ii] = 0;
	  }
	}

	RComplex<REALT> E_minus;
	RComplex<REALT> B_minus;
	RComplex<REALT> ctmp_0;
	RComplex<REALT> ctmp_1;
	RScalar<REALT> rtmp_0;
	RScalar<REALT> rtmp_1;

	for (int i = 0; i < Nc; ++i)
	{

	  /*# diag_L(i,0) = 1 - i*diag(E_z - B_z) */
	  /*#             = 1 - i*diag(F(3,2) - F(1,0)) */
	  ctmp_0 = f5.elem(site).elem().elem(i, i);
	  ctmp_0 -= f0.elem(site).elem().elem(i, i);
	  rtmp_0 = imag(ctmp_0);
	  tri[site].A.diag[0][i] += rtmp_0;

	  /*# diag_L(i+Nc,0) = 1 + i*diag(E_z - B_z) */
	  /*#                = 1 + i*diag(F(3,2) - F(1,0)) */
	  tri[site].A.diag[0][i + Nc] -= rtmp_0;

	  /*# diag_L(i,1) = 1 + i*diag(E_z + B_z) */
	  /*#             = 1 + i*diag(F(3,2) + F(1,0)) */
	  ctmp_1 = f5.elem(site).elem().elem(i, i);
	  ctmp_1 += f0.elem(site).elem().elem(i, i);
	  rtmp_1 = imag(ctmp_1);
	  tri[site].A.diag[1][i] -= rtmp_1;

	  /*# diag_L(i+Nc,1) = 1 - i*diag(E_z + B_z) */
	  /*#                = 1 - i*diag(F(3,2) + F(1,0)) */
	  tri[site].A.diag[1][i + Nc] += rtmp_1;
	}

	/*# Construct lower triangular portion */
	/*# Block diagonal terms */
	for (int i = 1; i < Nc; ++i)
	{
	  for (int j = 0; j < i; ++j)
	  {

	    int elem_ij = i * (i - 1) / 2 + j;
	    int elem_tmp = (i + Nc) * (i + Nc - 1) / 2 + j + Nc;

	    /*# L(i,j,0) = -i*(E_z - B_z)[i,j] */
	    /*#          = -i*(F(3,2) - F(1,0)) */
	    ctmp_0 = f0.elem(site).elem().elem(i, j);
	    ctmp_0 -= f5.elem(site).elem().elem(i, j);
	    tri[site].A.offd[0][elem_ij] = timesI(ctmp_0);

	    /*# L(i+Nc,j+Nc,0) = +i*(E_z - B_z)[i,j] */
	    /*#                = +i*(F(3,2) - F(1,0)) */
	    tri[site].A.offd[0][elem_tmp] = -tri[site].A.offd[0][elem_ij];

	    /*# L(i,j,1) = i*(E_z + B_z)[i,j] */
	    /*#          = i*(F(3,2) + F(1,0)) */
	    ctmp_1 = f5.elem(site).elem().elem(i, j);
	    ctmp_1 += f0.elem(site).elem().elem(i, j);
	    tri[site].A.offd[1][elem_ij] = timesI(ctmp_1);

	    /*# L(i+Nc,j+Nc,1) = -i*(E_z + B_z)[i,j] */
	    /*#                = -i*(F(3,2) + F(1,0)) */
	    tri[site].A.offd[1][elem_tmp] = -tri[site].A.offd[1][elem_ij];
	  }
	}

	/*# Off-diagonal */
	for (int i = 0; i < Nc; ++i)
	{
	  for (int j = 0; j < Nc; ++j)
	  {

	    // Flipped index
	    // by swapping i <-> j. In the past i would run slow
	    // and now j runs slow
	    int elem_ij = (i + Nc) * (i + Nc - 1) / 2 + j;

	    /*# i*E_- = (i*E_x + E_y) */
	    /*#       = (i*F(3,0) + F(3,1)) */
	    E_minus = timesI(f2.elem(site).elem().elem(i, j));
	    E_minus += f4.elem(site).elem().elem(i, j);

	    /*# i*B_- = (i*B_x + B_y) */
	    /*#       = (i*F(2,1) - F(2,0)) */
	    B_minus = timesI(f3.elem(site).elem().elem(i, j));
	    B_minus -= f1.elem(site).elem().elem(i, j);

	    /*# L(i+Nc,j,0) = -i*(E_- - B_-)  */
	    tri[site].A.offd[0][elem_ij] = B_minus - E_minus;

	    /*# L(i+Nc,j,1) = +i*(E_- + B_-)  */
	    tri[site].A.offd[1][elem_ij] = E_minus + B_minus;
	  }
	}
      } /* End Site loop */
#endif
    }	/* Function */
  }

  /* This now just sets up and dispatches... */
  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::makeClov(const multi1d<U>& f, const RealT& diag_mass)
  {
    START_CODE();

    if (Nd != 4)
    {
      QDPIO::cerr << __func__ << ": expecting Nd==4" << std::endl;
      QDP_abort(1);
    }

    if (Ns != 4)
    {
      QDPIO::cerr << __func__ << ": expecting Ns==4" << std::endl;
      QDP_abort(1);
    }

    U f0 = f[0] * getCloverCoeff(0, 1);
    U f1 = f[1] * getCloverCoeff(0, 2);
    U f2 = f[2] * getCloverCoeff(0, 3);
    U f3 = f[3] * getCloverCoeff(1, 2);
    U f4 = f[4] * getCloverCoeff(1, 3);
    U f5 = f[5] * getCloverCoeff(2, 3);

    const int nodeSites = QDP::Layout::sitesOnNode();
    QDPExpCloverEnv::QDPCloverMakeClovArg<U> arg = {diag_mass, f0, f1, f2, f3, f4, f5, tri};
    dispatch_to_threads(nodeSites, arg, QDPExpCloverEnv::makeClovSiteLoop<U>);

    END_CODE();
  }

  //! Invert
  /*!
   * Computes the inverse of the term on cb using Cholesky
   */
  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::choles(int cb)
  {
    START_CODE();

    // When you are doing the cholesky - also fill out the trace_log_diag piece)
    // chlclovms(tr_log_diag_, cb);
    // Switch to LDL^\dag inversion
    ldagdlinv(tr_M, cb);

    END_CODE();
  }

  //! Invert
  /*!
   * Computes the inverse of the term on cb using Cholesky
   *
   * \return logarithm of the determinant  
   */
  template <typename T, typename U, int N_exp>
  Double QDPExpCloverTermT<T, U, N_exp>::cholesDet(int cb) const
  {
#ifndef QDP_IS_QDPJIT
    return sum(tr_M, rb[cb]);
#else
    assert(!"ni");
    Double ret = 0.;
    return ret;
#endif
  }

  namespace QDPExpCloverEnv
  {

    template <typename U>
    struct LDagDLInvArgs {
      typedef typename WordType<U>::Type_t REALT;
      typedef OScalar<PScalar<PScalar<RScalar<REALT>>>> RealT;
      typedef OLattice<PScalar<PScalar<RScalar<REALT>>>> LatticeRealT;
      LatticeRealT& tr_Minv;
      ExpClovTriang<REALT>* tri;
      int cb;
    };

    template <typename U>
    inline void LDagDLInvSiteLoop(int lo, int hi, int myId, LDagDLInvArgs<U>* a)
    {
      typedef typename LDagDLInvArgs<U>::REALT REALT;
      typedef typename LDagDLInvArgs<U>::RealT RealT;
      typedef typename LDagDLInvArgs<U>::LatticeRealT LatticeRealT;

      LatticeRealT& tr_Minv = a->tr_Minv;
      ExpClovTriang<REALT>* tri = a->tri;
      int cb = a->cb;

      RScalar<REALT> zip = 0;
      int N = 2 * Nc;

      // Loop through the sites.
      for (int ssite = lo; ssite < hi; ++ssite)
      {

	int site = rb[cb].siteTable()[ssite];
	for (int block = 0; block < 2; block++)
	{
	  for (int j = 0; j < 6; ++j)
	  {
	    auto tmp = tri[site].q[block][j];
	    tri[site].q[block][j] = tri[site].qinv[block][j];
	    tri[site].qinv[block][j] = tmp;
	  }
	}

      } /* End Site Loop */
    }	/* End Function */

    template <typename U>
    struct TracePowersArgs {
      typedef typename WordType<U>::Type_t REALT;
      typedef OScalar<PScalar<PScalar<RScalar<REALT>>>> RealT;
      typedef OLattice<PScalar<PScalar<RScalar<REALT>>>> LatticeRealT;

      ExpClovTriang<REALT>* tri;
      int cb;
    };

    template <typename U, int N_exp>
    inline void TracePowerSiteLoop(int lo, int hi, int myId, TracePowersArgs<U>* a)
    {
      typedef typename LDagDLInvArgs<U>::REALT REALT;
      typedef typename LDagDLInvArgs<U>::RealT RealT;
      typedef typename LDagDLInvArgs<U>::LatticeRealT LatticeRealT;

      ExpClovTriang<REALT>* tri = a->tri;
      int cb = a->cb;

      RScalar<REALT> zip = 0;

      // Loop through the sites.
      for (int ssite = lo; ssite < hi; ++ssite)
      {

	int site = rb[cb].siteTable()[ssite];

	// Compute the q-s for the block
	REALT tab[2][N_exp + 1][6] = {0};
	REALT itab[2][N_exp + 1][6] = {0};

	for (int block = 0; block < 2; ++block)
	{
	  constexpr int upper = N_exp + 1 < 6 ? N_exp + 1 : 6;
	  REALT ifact = 1;
	  for (int i = 0; i < upper; ++i)
	  {
	    tab[block][i][i] = (REALT)1;
	    itab[block][i][i] = ifact;
	    ifact = -ifact;
	  }
	}

	Traces<REALT, 0> tr0(tri[site]);
	Traces<REALT, 1> tr1(tri[site]);

	multi1d<RScalar<REALT>> trace0(6);
	multi1d<RScalar<REALT>> trace1(6);
	tr0.traces(trace0);
	tr1.traces(trace1);

	/*
          trace[1][0] = toDouble(tr1.trace());
          trace[1][1] = toDouble(tr1.trace2());
          trace[1][2] = toDouble(tr1.trace3());
          trace[1][3] = toDouble(tr1.trace4());
          trace[1][4] = toDouble(tr1.trace5());
          trace[1][5] = toDouble(tr1.trace6());
          */
	if (N_exp + 1 > 6)
	{
	  REALT trace[2][6];
	  for (int i = 0; i < 6; ++i)
	    trace[0][i] = toDouble(trace0[i]);
	  for (int i = 0; i < 6; ++i)
	    trace[1][i] = toDouble(trace1[i]);

	  for (int block = 0; block < 2; ++block)
	  {
	    REALT p[5];

	    p[4] = ((REALT)1 / (REALT)2) * trace[block][1]; // (1/2) Tr A^2

	    p[3] = ((REALT)1 / (REALT)3) * trace[block][2]; // (1/3) Tr A^3
	    p[2] = ((REALT)1 / (REALT)4) * trace[block][3] -
		   ((REALT)1 / (REALT)8) * trace[block][1] *
		     trace[block][1]; // (1/4) Tr A^4 - (1/8) (Tr A^2)^2

	    p[1] = ((REALT)1 / (REALT)5) * trace[block][4] -
		   ((REALT)1 / (REALT)6) * trace[block][2] *
		     trace[block][1];			   // (1/5) Tr A^5 - (1/6) Tr A^3 Tr A^2

	    p[0] = ((REALT)1 / (REALT)6) * trace[block][5] // (1/6) Tr A^6 - (1/8) Tr A^4 Tr A^2
		   - ((REALT)1 / (REALT)8) * trace[block][3] *
		       trace[block][1]			   //     - (1/18) [ Tr A^3 ]^2
		   - ((REALT)1 / (REALT)18) * trace[block][2] *
		       trace[block][2]			   //     + (1/48) [ Tr A^2 ]^3
		   + ((REALT)1 / (REALT)48) * trace[block][1] * trace[block][1] * trace[block][1];

	    // Row 6
	    for (int i = 0; i < 5; ++i)
	    {
	      tab[block][6][i] = p[i];
	    }

	    // Row 7+

	    for (int row = 7; row <= N_exp; ++row)
	    {
	      for (int i = 0; i < 5; i++)
	      {
		for (int j = 0; j < 6; j++)
		{
		  tab[block][row][j] += p[i] * tab[block][row - 6 + i][j];
		}
	      }
	    }
	  } // Blocks
	}   // N_exp + 1 >

	// Sum into the q
	for (int block = 0; block < 2; ++block)
	{

	  // Row 0
	  for (int i = 0; i < 6; i++)
	  {
	    tri[site].q[block][i] = RScalar<REALT>(tab[block][0][i]);
	    tri[site].qinv[block][i] = RScalar<REALT>(tab[block][0][i]);
	  }

	  unsigned long fact = 1;
	  for (unsigned int row = 1; row <= N_exp; ++row)
	  {
	    fact *= (unsigned long)row;
	    REALT sign = (row % 2 == 0) ? (REALT)1 : (REALT)(-1);
	    for (int i = 0; i < 6; i++)
	    {
	      tri[site].q[block][i] += RScalar<REALT>(tab[block][row][i] / (REALT)(fact));
	      tri[site].qinv[block][i] += RScalar<REALT>(sign * tab[block][row][i] / (REALT)(fact));
	    }
	  }
	}

	// Assemble te exponential from the q-s and powers of A.
	// siteExponentiate(tri[site]);
      } /* End Site Loop */
    }	/* End Function */

  }	/* End Namespace */

  /*! An LDL^\dag decomposition and inversion? */
  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::ldagdlinv(LatticeREAL& tr_Minv, int cb)
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    if (2 * Nc < 3)
    {
      QDPIO::cerr << __func__ << ": Matrix is too small" << std::endl;
      QDP_abort(1);
    }

    // Zero trace log
    tr_Minv[rb[cb]] = zero;

    QDPExpCloverEnv::LDagDLInvArgs<U> a = {tr_Minv, tri, cb};
    int num_site_table = rb[cb].numSiteTable();
    dispatch_to_threads(num_site_table, a, QDPExpCloverEnv::LDagDLInvSiteLoop<U>);

    END_CODE();
#endif
  }

  /* This now just sets up and dispatches... */
  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::tracePowers()
  {
    START_CODE();

    if (Nd != 4)
    {
      QDPIO::cerr << __func__ << ": expecting Nd==4" << std::endl;
      QDP_abort(1);
    }

    if (Ns != 4)
    {
      QDPIO::cerr << __func__ << ": expecting Ns==4" << std::endl;
      QDP_abort(1);
    }

    for (int cb = 0; cb < 2; ++cb)
    {
      const int num_sites = rb[cb].siteTable().size();
      QDPExpCloverEnv::TracePowersArgs<U> arg = {tri, cb};
      dispatch_to_threads(num_sites, arg, QDPExpCloverEnv::TracePowerSiteLoop<U, N_exp>);
    }
    END_CODE();
  }

  /*! CHLCLOVMS - Cholesky decompose the clover mass term and uses it to
   *              compute  lower(A^-1) = lower((L.L^dag)^-1)
   *              Adapted from Golub and Van Loan, Matrix Computations, 2nd, Sec 4.2.4
   *
   * Arguments:
   *
   * \param DetP         flag whether to compute determinant (Read)
   * \param logdet       logarithm of the determinant        (Write)
   * \param cb           checkerboard of work                (Read)
   */

  /**
   * Apply a dslash
   *
   * Performs the operation
   *
   *  chi <-   (L + D + L^dag) . psi
   *
   * where
   *   L       is a lower triangular matrix
   *   D       is the real diagonal. (stored together in type TRIANG)
   *
   * Arguments:
   * \param chi     result                                      (Write)
   * \param psi     source                                      (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
   */
  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::applyPowerSite(T& chi, const T& psi, enum PlusMinus isign,
						      int site, int power) const
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    if (Ns != 4)
    {
      QDPIO::cerr << __func__ << ": CloverTerm::applySite requires Ns==4" << std::endl;
      QDP_abort(1);
    }

    RComplex<REALT>* cchi = (RComplex<REALT>*)&(chi.elem(site).elem(0).elem(0));
    const RComplex<REALT>* const ppsi =
      (const RComplex<REALT>* const)&(psi.elem(site).elem(0).elem(0));

    siteApplicationPower<REALT>(cchi, tri, ppsi, power);

    END_CODE();
#endif
  }

  /**
   * Apply a dslash
   *
   * Performs the operation
   *
   *  chi <-   (L + D + L^dag) . psi
   *
   * where
   *   L       is a lower triangular matrix
   *   D       is the real diagonal. (stored together in type TRIANG)
   *
   * Arguments:
   * \param chi     result                                      (Write)
   * \param psi     source                                      (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
   */
  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::applySite(T& chi, const T& psi, enum PlusMinus isign,
						 int site) const
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    if (Ns != 4)
    {
      QDPIO::cerr << __func__ << ": CloverTerm::applySite requires Ns==4" << std::endl;
      QDP_abort(1);
    }

    RComplex<REALT>* cchi = (RComplex<REALT>*)&(chi.elem(site).elem(0).elem(0));
    const RComplex<REALT>* const ppsi =
      (const RComplex<REALT>* const)&(psi.elem(site).elem(0).elem(0));
    siteApplicationExp(cchi, tri[site], ppsi);
    END_CODE();
#endif
  }

  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::applySiteInv(T& chi, const T& psi, enum PlusMinus isign,
						    int site) const
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    if (Ns != 4)
    {
      QDPIO::cerr << __func__ << ": CloverTerm::applySite requires Ns==4" << std::endl;
      QDP_abort(1);
    }

    RComplex<REALT>* cchi = (RComplex<REALT>*)&(chi.elem(site).elem(0).elem(0));
    const RComplex<REALT>* const ppsi =
      (const RComplex<REALT>* const)&(psi.elem(site).elem(0).elem(0));
    siteApplicationExp<REALT, 1>(cchi, tri[site], ppsi);
    END_CODE();
#endif
  }

  //! Returns the appropriate clover coefficient for indices mu and nu
  template <typename T, typename U, int N_exp>
  Real QDPExpCloverTermT<T, U, N_exp>::getCloverCoeff(int mu, int nu) const
  {
    START_CODE();

    if (param.anisoParam.anisoP)
    {
      if (mu == param.anisoParam.t_dir || nu == param.anisoParam.t_dir)
      {
	return param.clovCoeffT;
      }
      else
      {
	// Otherwise return the spatial coeff
	return param.clovCoeffR;
      }
    }
    else
    {
      // If there is no anisotropy just return the spatial one, it will
      // be the same as the temporal one
      return param.clovCoeffR;
    }

    END_CODE();
  }

  namespace QDPExpCloverEnv
  {
    template <typename T>
    struct ApplyPowerArgs {
      typedef typename WordType<T>::Type_t REALT;
      T& chi;
      const T& psi;
      const ExpClovTriang<REALT>* tri;
      int cb;
      int power = 1;
    };

    template <typename T>
    void applySitePowerLoop(int lo, int hi, int MyId, ApplyPowerArgs<T>* arg)
    {
#ifndef QDP_IS_QDPJIT
      // This is essentially the body of the previous "Apply"
      // but now the args are handed in through user arg struct...

      START_CODE();

      typedef typename WordType<T>::Type_t REALT;
      // Unwrap the args...
      T& chi = arg->chi;
      const T& psi = arg->psi;
      const ExpClovTriang<REALT>* tri = arg->tri;
      int cb = arg->cb;
      int power = arg->power;
      const int n = 2 * Nc;

      for (int ssite = lo; ssite < hi; ++ssite)
      {

	int site = rb[cb].siteTable()[ssite];

	RComplex<REALT>* cchi = (RComplex<REALT>*)&(chi.elem(site).elem(0).elem(0));

	const RComplex<REALT>* const ppsi =
	  (const RComplex<REALT>* const)&(psi.elem(site).elem(0).elem(0));

	siteApplicationPower<REALT>(cchi, tri[site], ppsi, power);
      }
      END_CODE();
#endif
    } // Function

    template <typename T>
    struct ApplyArgs {
      typedef typename WordType<T>::Type_t REALT;
      T& chi;
      const T& psi;
      const ExpClovTriang<REALT>* tri;
      int cb;
    };

    template <typename T, int inv = 0>
    void applySiteLoop(int lo, int hi, int MyId, ApplyArgs<T>* arg)
    {
#ifndef QDP_IS_QDPJIT
      // This is essentially the body of the previous "Apply"
      // but now the args are handed in through user arg struct...

      START_CODE();

      typedef typename WordType<T>::Type_t REALT;
      // Unwrap the args...
      T& chi = arg->chi;
      const T& psi = arg->psi;
      const ExpClovTriang<REALT>* tri = arg->tri;
      int cb = arg->cb;
      const int n = 2 * Nc;

      for (int ssite = lo; ssite < hi; ++ssite)
      {

	int site = rb[cb].siteTable()[ssite];
	RComplex<REALT>* cchi = (RComplex<REALT>*)&(chi.elem(site).elem(0).elem(0));
	const RComplex<REALT>* const ppsi =
	  (const RComplex<REALT>* const)&(psi.elem(site).elem(0).elem(0));

	siteApplicationExp<REALT, inv>(cchi, tri[site], ppsi);
      }
      END_CODE();
#endif
    } // Function
  }   // Namespace

  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::applyRef(T& chi, const T& psi, enum PlusMinus isign,
						int N) const
  {
    chi = psi;
    T tmp;
    for (int j = N; j > 1; --j)
    {
      for (int cb = 0; cb < 2; ++cb)
	(*this).applyPower(tmp, chi, isign, cb, 1); // M^1 chi
      tmp /= Real(j);				    // (1/N) M chi
      chi = psi + tmp;				    // ( psi + (1/N)M chi)
    }
    for (int cb = 0; cb < 2; ++cb)
      (*this).applyPower(tmp, chi, isign, cb, 1); // M chi

    chi = psi + tmp;				  // psi + M chi
    (*this).getFermBC().modifyF(chi);
  }
  /**
   * Apply a dslash
   *
   * Performs the operation
   *
   *  chi <-   (L + D + L^dag) . psi
   *
   * where
   *   L       is a lower triangular matrix
   *   D       is the real diagonal. (stored together in type TRIANG)
   *
   * Arguments:
   * \param chi     result                                      (Write)
   * \param psi     source                                      (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of OUTPUT std::vector               (Read) 
   */

  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::applyPower(T& chi, const T& psi, enum PlusMinus isign,
						  int cb, int power) const
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    if (Ns != 4)
    {
      QDPIO::cerr << __func__ << ": CloverTerm::apply requires Ns==4" << std::endl;
      QDP_abort(1);
    }

    QDPExpCloverEnv::ApplyPowerArgs<T> arg = {chi, psi, tri, cb, power};
    int num_sites = rb[cb].siteTable().size();

    // The dispatch function is at the end of the file
    // ought to work for non-threaded targets too...
    dispatch_to_threads(num_sites, arg, QDPExpCloverEnv::applySitePowerLoop<T>);
    (*this).getFermBC().modifyF(chi, QDP::rb[cb]);

    END_CODE();
#endif
  }

  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::apply(T& chi, const T& psi, enum PlusMinus isign,
					     int cb) const
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    if (Ns != 4)
    {
      QDPIO::cerr << __func__ << ": CloverTerm::apply requires Ns==4" << std::endl;
      QDP_abort(1);
    }

    QDPExpCloverEnv::ApplyArgs<T> arg = {chi, psi, tri, cb};
    int num_sites = rb[cb].siteTable().size();

    // The dispatch function is at the end of the file
    // ought to work for non-threaded targets too...
    dispatch_to_threads(num_sites, arg, QDPExpCloverEnv::applySiteLoop<T, 0>);
    (*this).getFermBC().modifyF(chi, QDP::rb[cb]);

    END_CODE();
#endif
  }

  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::applyInv(T& chi, const T& psi, enum PlusMinus isign,
						int cb) const
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    if (Ns != 4)
    {
      QDPIO::cerr << __func__ << ": CloverTerm::apply requires Ns==4" << std::endl;
      QDP_abort(1);
    }

    QDPExpCloverEnv::ApplyArgs<T> arg = {chi, psi, tri, cb};
    int num_sites = rb[cb].siteTable().size();

    // The dispatch function is at the end of the file
    // ought to work for non-threaded targets too...
    dispatch_to_threads(num_sites, arg, QDPExpCloverEnv::applySiteLoop<T, 1>);
    (*this).getFermBC().modifyF(chi, QDP::rb[cb]);

    END_CODE();
#endif
  }

  //! TRIACNTR
  /*! 
   * \ingroup linop
   *
   *  Calculates
   *     Tr_D ( Gamma_mat L )
   *
   * This routine is specific to Wilson fermions!
   * 
   *  the trace over the Dirac indices for one of the 16 Gamma matrices
   *  and a hermitian color x spin matrix A, stored as a block diagonal
   *  complex lower triangular matrix L and a real diagonal diag_L.

   *  Here 0 <= mat <= 15 and
   *  if mat = mat_1 + mat_2 * 2 + mat_3 * 4 + mat_4 * 8
   *
   *  Gamma(mat) = gamma(1)^(mat_1) * gamma(2)^(mat_2) * gamma(3)^(mat_3)
   *             * gamma(4)^(mat_4)
   *
   *  Further, in basis for the Gamma matrices used, A is of the form
   *
   *      | A_0 |  0  |
   *  A = | --------- |
   *      |  0  | A_1 |
   *
   *
   * Arguments:
   *
   *  \param B         the resulting SU(N) color matrix	  (Write) 
   *  \param clov      clover term                        (Read) 
   *  \param mat       label of the Gamma matrix          (Read)
   */

  namespace QDPExpCloverEnv
  {
    template <typename U>
    struct TriaCntrArgs {
      typedef typename WordType<U>::Type_t REALT;

      U& B;
      const ExpClovTriang<REALT>* tri;
      int mat;
      int cb;
    };

    template <typename U>
    inline void triaCntrSiteLoop(int lo, int hi, int myId, TriaCntrArgs<U>* a)
    {
      typedef typename WordType<U>::Type_t REALT;
      U& B = a->B;
      const ExpClovTriang<REALT>* tri = a->tri;
      int mat = a->mat;
      int cb = a->cb;

      for (int ssite = lo; ssite < hi; ++ssite)
      {

	int site = rb[cb].siteTable()[ssite];

	switch (mat)
	{

	case 0:
	  /*# gamma(   0)   1  0  0  0            # ( 0000 )  --> 0 */
	  /*#               0  1  0  0 */
	  /*#               0  0  1  0 */
	  /*#               0  0  0  1 */
	  /*# From diagonal part */
	  {
	    RComplex<REALT> lctmp0;
	    RScalar<REALT> lr_zero0;
	    RScalar<REALT> lrtmp0;

	    lr_zero0 = 0;

	    for (int i0 = 0; i0 < Nc; ++i0)
	    {

	      lrtmp0 = tri[site].A.diag[0][i0];
	      lrtmp0 += tri[site].A.diag[0][i0 + Nc];
	      lrtmp0 += tri[site].A.diag[1][i0];
	      lrtmp0 += tri[site].A.diag[1][i0 + Nc];
	      B.elem(site).elem().elem(i0, i0) = cmplx(lrtmp0, lr_zero0);
	    }

	    /*# From lower triangular portion */
	    int elem_ij0 = 0;
	    for (int i0 = 1; i0 < Nc; ++i0)
	    {

	      int elem_ijb0 = (i0 + Nc) * (i0 + Nc - 1) / 2 + Nc;

	      for (int j0 = 0; j0 < i0; ++j0)
	      {

		lctmp0 = tri[site].A.offd[0][elem_ij0];
		lctmp0 += tri[site].A.offd[0][elem_ijb0];
		lctmp0 += tri[site].A.offd[1][elem_ij0];
		lctmp0 += tri[site].A.offd[1][elem_ijb0];

		B.elem(site).elem().elem(j0, i0) = lctmp0;
		B.elem(site).elem().elem(i0, j0) = adj(lctmp0);

		elem_ij0++;
		elem_ijb0++;
	      }
	    }
	  }

	  break;

	case 3:
	  /*# gamma(  12)  -i  0  0  0            # ( 0011 )  --> 3 */
	  /*#               0  i  0  0 */
	  /*#               0  0 -i  0 */
	  /*#               0  0  0  i */
	  /*# From diagonal part */

	  {
	    RComplex<REALT> lctmp3;
	    RScalar<REALT> lr_zero3;
	    RScalar<REALT> lrtmp3;

	    lr_zero3 = 0;

	    for (int i3 = 0; i3 < Nc; ++i3)
	    {

	      lrtmp3 = tri[site].A.diag[0][i3 + Nc];
	      lrtmp3 -= tri[site].A.diag[0][i3];
	      lrtmp3 -= tri[site].A.diag[1][i3];
	      lrtmp3 += tri[site].A.diag[1][i3 + Nc];
	      B.elem(site).elem().elem(i3, i3) = cmplx(lr_zero3, lrtmp3);
	    }

	    /*# From lower triangular portion */
	    int elem_ij3 = 0;
	    for (int i3 = 1; i3 < Nc; ++i3)
	    {

	      int elem_ijb3 = (i3 + Nc) * (i3 + Nc - 1) / 2 + Nc;

	      for (int j3 = 0; j3 < i3; ++j3)
	      {

		lctmp3 = tri[site].A.offd[0][elem_ijb3];
		lctmp3 -= tri[site].A.offd[0][elem_ij3];
		lctmp3 -= tri[site].A.offd[1][elem_ij3];
		lctmp3 += tri[site].A.offd[1][elem_ijb3];

		B.elem(site).elem().elem(j3, i3) = timesI(adj(lctmp3));
		B.elem(site).elem().elem(i3, j3) = timesI(lctmp3);

		elem_ij3++;
		elem_ijb3++;
	      }
	    }
	  }
	  break;

	case 5:
	  /*# gamma(  13)   0 -1  0  0            # ( 0101 )  --> 5 */
	  /*#               1  0  0  0 */
	  /*#               0  0  0 -1 */
	  /*#               0  0  1  0 */

	  {

	    RComplex<REALT> lctmp5;
	    RScalar<REALT> lrtmp5;

	    for (int i5 = 0; i5 < Nc; ++i5)
	    {

	      int elem_ij5 = (i5 + Nc) * (i5 + Nc - 1) / 2;

	      for (int j5 = 0; j5 < Nc; ++j5)
	      {

		int elem_ji5 = (j5 + Nc) * (j5 + Nc - 1) / 2 + i5;

		lctmp5 = adj(tri[site].A.offd[0][elem_ji5]);
		lctmp5 -= tri[site].A.offd[0][elem_ij5];
		lctmp5 += adj(tri[site].A.offd[1][elem_ji5]);
		lctmp5 -= tri[site].A.offd[1][elem_ij5];

		B.elem(site).elem().elem(i5, j5) = lctmp5;

		elem_ij5++;
	      }
	    }
	  }
	  break;

	case 6:
	  /*# gamma(  23)   0 -i  0  0            # ( 0110 )  --> 6 */
	  /*#              -i  0  0  0 */
	  /*#               0  0  0 -i */
	  /*#               0  0 -i  0 */

	  {
	    RComplex<REALT> lctmp6;
	    RScalar<REALT> lrtmp6;

	    for (int i6 = 0; i6 < Nc; ++i6)
	    {

	      int elem_ij6 = (i6 + Nc) * (i6 + Nc - 1) / 2;

	      for (int j6 = 0; j6 < Nc; ++j6)
	      {

		int elem_ji6 = (j6 + Nc) * (j6 + Nc - 1) / 2 + i6;

		lctmp6 = adj(tri[site].A.offd[0][elem_ji6]);
		lctmp6 += tri[site].A.offd[0][elem_ij6];
		lctmp6 += adj(tri[site].A.offd[1][elem_ji6]);
		lctmp6 += tri[site].A.offd[1][elem_ij6];

		B.elem(site).elem().elem(i6, j6) = timesMinusI(lctmp6);

		elem_ij6++;
	      }
	    }
	  }
	  break;

	case 9:
	  /*# gamma(  14)   0  i  0  0            # ( 1001 )  --> 9 */
	  /*#               i  0  0  0 */
	  /*#               0  0  0 -i */
	  /*#               0  0 -i  0 */

	  {
	    RComplex<REALT> lctmp9;
	    RScalar<REALT> lrtmp9;

	    for (int i9 = 0; i9 < Nc; ++i9)
	    {

	      int elem_ij9 = (i9 + Nc) * (i9 + Nc - 1) / 2;

	      for (int j9 = 0; j9 < Nc; ++j9)
	      {

		int elem_ji9 = (j9 + Nc) * (j9 + Nc - 1) / 2 + i9;

		lctmp9 = adj(tri[site].A.offd[0][elem_ji9]);
		lctmp9 += tri[site].A.offd[0][elem_ij9];
		lctmp9 -= adj(tri[site].A.offd[1][elem_ji9]);
		lctmp9 -= tri[site].A.offd[1][elem_ij9];

		B.elem(site).elem().elem(i9, j9) = timesI(lctmp9);

		elem_ij9++;
	      }
	    }
	  }
	  break;

	case 10:
	  /*# gamma(  24)   0 -1  0  0            # ( 1010 )  --> 10 */
	  /*#               1  0  0  0 */
	  /*#               0  0  0  1 */
	  /*#               0  0 -1  0 */
	  {
	    RComplex<REALT> lctmp10;
	    RScalar<REALT> lrtmp10;

	    for (int i10 = 0; i10 < Nc; ++i10)
	    {

	      int elem_ij10 = (i10 + Nc) * (i10 + Nc - 1) / 2;

	      for (int j10 = 0; j10 < Nc; ++j10)
	      {

		int elem_ji10 = (j10 + Nc) * (j10 + Nc - 1) / 2 + i10;

		lctmp10 = adj(tri[site].A.offd[0][elem_ji10]);
		lctmp10 -= tri[site].A.offd[0][elem_ij10];
		lctmp10 -= adj(tri[site].A.offd[1][elem_ji10]);
		lctmp10 += tri[site].A.offd[1][elem_ij10];

		B.elem(site).elem().elem(i10, j10) = lctmp10;

		elem_ij10++;
	      }
	    }
	  }
	  break;

	case 12:
	  /*# gamma(  34)   i  0  0  0            # ( 1100 )  --> 12 */
	  /*#               0 -i  0  0 */
	  /*#               0  0 -i  0 */
	  /*#               0  0  0  i */
	  /*# From diagonal part */
	  {

	    RComplex<REALT> lctmp12;
	    RScalar<REALT> lr_zero12;
	    RScalar<REALT> lrtmp12;

	    lr_zero12 = 0;

	    for (int i12 = 0; i12 < Nc; ++i12)
	    {

	      lrtmp12 = tri[site].A.diag[0][i12];
	      lrtmp12 -= tri[site].A.diag[0][i12 + Nc];
	      lrtmp12 -= tri[site].A.diag[1][i12];
	      lrtmp12 += tri[site].A.diag[1][i12 + Nc];
	      B.elem(site).elem().elem(i12, i12) = cmplx(lr_zero12, lrtmp12);
	    }

	    /*# From lower triangular portion */
	    int elem_ij12 = 0;
	    for (int i12 = 1; i12 < Nc; ++i12)
	    {

	      int elem_ijb12 = (i12 + Nc) * (i12 + Nc - 1) / 2 + Nc;

	      for (int j12 = 0; j12 < i12; ++j12)
	      {

		lctmp12 = tri[site].A.offd[0][elem_ij12];
		lctmp12 -= tri[site].A.offd[0][elem_ijb12];
		lctmp12 -= tri[site].A.offd[1][elem_ij12];
		lctmp12 += tri[site].A.offd[1][elem_ijb12];

		B.elem(site).elem().elem(i12, j12) = timesI(lctmp12);
		B.elem(site).elem().elem(j12, i12) = timesI(adj(lctmp12));

		elem_ij12++;
		elem_ijb12++;
	      }
	    }
	  }
	  break;

	default: QDPIO::cout << __func__ << ": invalid Gamma matrix int" << std::endl; QDP_abort(1);
	}

      } // END Site Loop
    }	// End Function
  }	// End Namespace

  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::triacntr(U& B, int mat, int cb) const
  {
#ifndef QDP_IS_QDPJIT
    START_CODE();

    B = zero;

    if (mat < 0 || mat > 15)
    {
      QDPIO::cerr << __func__ << ": Gamma out of range: mat = " << mat << std::endl;
      QDP_abort(1);
    }

    QDPExpCloverEnv::TriaCntrArgs<U> a = {B, tri, mat, cb};
    dispatch_to_threads(rb[cb].numSiteTable(), a, QDPExpCloverEnv::triaCntrSiteLoop<U>);

    END_CODE();
#endif
  }

  namespace QDPExpCloverEnv
  {
    template <typename R>
    struct QUDAPackArgs {
      int cb;
      multi1d<QUDAPackedClovSite<R>>& quda_array;
      const ExpClovTriang<R>* tri;
    };

    template <typename R>
    void qudaPackSiteLoop(int lo, int hi, int myId, QUDAPackArgs<R>* a)
    {
#if 0
      int cb = a->cb;
      int Ns2 = Ns / 2;

      multi1d<QUDAPackedClovSite<R>>& quda_array = a->quda_array;
      const ExpClovTriang<R>* tri = a->tri;

      const int idtab[15] = {0, 1, 3, 6, 10, 2, 4, 7, 11, 5, 8, 12, 9, 13, 14};

      for (int ssite = lo; ssite < hi; ++ssite) {
      
	      int site = rb[cb].siteTable()[ssite];
	      // First Chiral Block
	      for (int i = 0; i < 6; i++) {
      	  quda_array[site].diag1[i] = tri[site].Exp[0].diag[0][i].elem();
	      }

      	int target_index = 0;

	      for (int col = 0; col < Nc * Ns2 - 1; col++) {
	        for (int row = col + 1; row < Nc * Ns2; row++) {

    	      int source_index = row * (row - 1) / 2 + col;

	          quda_array[site].offDiag1[target_index][0] = tri[site].Exp[0].offd[0][source_index].real();
	          quda_array[site].offDiag1[target_index][1] = tri[site].Exp[0].offd[0][source_index].imag();
	          target_index++;
	        }
	      }


	      // Second Chiral Block
	      for (int i = 0; i < 6; i++) {
	        quda_array[site].diag2[i] = tri[site].Exp[0].diag[1][i].elem();
	      }

      	target_index = 0;
	      for (int col = 0; col < Nc * Ns2 - 1; col++) {
	        for (int row = col + 1; row < Nc * Ns2; row++) {

      	    int source_index = row * (row - 1) / 2 + col;
	          quda_array[site].offDiag2[target_index][0] = tri[site].Exp[0].offd[1][source_index].real();
	          quda_array[site].offDiag2[target_index][1] = tri[site].Exp[0].offd[1][source_index].imag();
	          target_index++;
	        }
	      }
      }
#endif
    }
  }

  template <typename T, typename U, int N_exp>
  void QDPExpCloverTermT<T, U, N_exp>::packForQUDA(
    multi1d<QUDAPackedClovSite<typename WordType<T>::Type_t>>& quda_array, int cb) const
  {
    typedef typename WordType<T>::Type_t REALT;
    int num_sites = rb[cb].siteTable().size();

    QDPExpCloverEnv::QUDAPackArgs<REALT> args = {cb, quda_array, tri};
    dispatch_to_threads(num_sites, args, QDPExpCloverEnv::qudaPackSiteLoop<REALT>);
  }

  template <int N = N_exp_default>
  using QDPExpCloverTerm = QDPExpCloverTermT<LatticeFermion, LatticeColorMatrix, N>;

  template <int N = N_exp_default>
  using QDPExpCloverTermF = QDPExpCloverTermT<LatticeFermionF, LatticeColorMatrixF, N>;

  template <int N = N_exp_default>
  using QDPExpCloverTermD = QDPExpCloverTermT<LatticeFermionD, LatticeColorMatrixD, N>;
} // End Namespace Chroma

#endif
