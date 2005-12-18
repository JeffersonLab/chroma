// $Id: clover_term_qdp_w.cc,v 2.1 2005-12-18 23:53:26 edwards Exp $
/*! \file
 *  \brief Clover term linear operator
 *
 *  This particular implementation is specific to a scalar-like
 *  architecture
 */

#include "chromabase.h"
#include "actions/ferm/linop/clover_term_qdp_w.h"


namespace Chroma 
{ 

  typedef PColorMatrix < RComplex <REAL>, Nc > PrimitiveSU3Matrix;
  typedef RComplex<REAL> PrimitiveLDUFerm[2][2*Nc];



  //! Creation routine
  void QDPCloverTerm::create(const multi1d<LatticeColorMatrix>& u_, 	
			     const CloverFermActParams& param_)
  {
    u = u_;
    param = param_;

    tri.resize(QDP::Layout::sitesOnNode());  // hold local lattice

    /* Calculate F(mu,nu) */
    multi1d<LatticeColorMatrix> f;
    mesField(f, u);
    makeClov(f);
  }


  /*
   * MAKCLOV 
   *
   *  In this routine, MAKCLOV calculates

   *    1 - (1/4)*sigma(mu,nu) F(mu,nu)

   *  using F from mesfield

   *    F(mu,nu) =  (1/4) sum_p (1/2) [ U_p(x) - U^dag_p(x) ]

   *  using basis of SPPROD and stores in a lower triangular matrix
   *  (no diagonal) plus real diagonal

   *  where
   *    U_1 = u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)
   *    U_2 = u(x,nu)*u_dag(x-mu+nu,mu)*u_dag(x-mu,nu)*u(x-mu,mu)
   *    U_3 = u_dag(x-mu,mu)*u_dag(x-mu-nu,nu)*u(x-mu-nu,mu)*u(x-nu,nu)
   *    U_4 = u_dag(x-nu,nu)*u(x-nu,mu)*u(x-nu+mu,nu)*u_dag(x,mu)

   *  and

   *         | sigF(1)   sigF(3)     0         0     |
   *  sigF = | sigF(5)  -sigF(1)     0         0     |
   *         |   0         0      -sigF(0)  -sigF(2) |
   *         |   0         0      -sigF(4)   sigF(0) |
   *  where
   *    sigF(i) is a color matrix

   *  sigF(0) = i*(ClovT*E_z + ClovR*B_z)
   *          = i*(ClovT*F(3,2) + ClovR*F(1,0))
   *  sigF(1) = i*(ClovT*E_z - ClovR*B_z)
   *          = i*(ClovT*F(3,2) - ClovR*F(1,0))
   *  sigF(2) = i*(E_+ + B_+)
   *  sigF(3) = i*(E_+ - B_+)
   *  sigF(4) = i*(E_- + B_-)
   *  sigF(5) = i*(E_- - B_-)
   *  i*E_+ = (i*ClovT*E_x - ClovT*E_y)
   *        = (i*ClovT*F(3,0) - ClovT*F(3,1))
   *  i*E_- = (i*ClovT*E_x + ClovT*E_y)
   *        = (i*ClovT*F(3,0) + ClovT*F(3,1))
   *  i*B_+ = (i*ClovR*B_x - ClovR*B_y)
   *        = (i*ClovR*F(2,1) + ClovR*F(2,0))
   *  i*B_- = (i*ClovR*B_x + ClovR*B_y)
   *        = (i*ClovR*F(2,1) - ClovR*F(2,0))

   *  NOTE: I am using  i*F  of the usual F defined by UKQCD, Heatlie et.al.

   *  NOTE: the above definitions assume that the time direction, t_dir,
   *        is 3. In general F(k,j) is multiplied with ClovT if either
   *        k=t_dir or j=t_dir, and with ClovR otherwise.

   *+++
   *  Here are some notes on the origin of this routine. NOTE, ClovCoeff or u0
   *  are not actually used in MAKCLOV.
   *
   *  The clover mass term is suppose to act on a vector like
   *
   *  chi = (1 - (ClovCoeff/u0^3) * kappa/4 * sum_mu sum_nu F(mu,nu)*sigma(mu,nu)) * psi

   *  Definitions used here (NOTE: no "i")
   *   sigma(mu,nu) = gamma(mu)*gamma(nu) - gamma(nu)*gamma(mu)
   *                = 2*gamma(mu)*gamma(nu)   for mu != nu
   *       
   *   chi = sum_mu sum_nu F(mu,nu)*gamma(mu)*gamma(nu)*psi   for mu < nu
   *       = (1/2) * sum_mu sum_nu F(mu,nu)*gamma(mu)*gamma(nu)*psi   for mu != nu
   *       = (1/4) * sum_mu sum_nu F(mu,nu)*sigma(mu,nu)*psi
   *
   *
   * chi = (1 - (ClovCoeff/u0^3) * kappa/4 * sum_mu sum_nu F(mu,nu)*sigma(mu,nu)) * psi
   *     = psi - (ClovCoeff/u0^3) * kappa * chi
   *     == psi - kappa * chi
   *
   *  We have absorbed ClovCoeff/u0^3 into kappa. A u0 was previously absorbed into kappa
   *  for compatibility to ancient conventions. 
   *---

   * Arguments:
   *  \param f         field strength tensor F(cb,mu,nu)        (Read)
   */
  void QDPCloverTerm::makeClov(const multi1d<LatticeColorMatrix>& f)
  {
    LatticeColorMatrix f0;
    LatticeColorMatrix f1;
    LatticeColorMatrix f2;
    LatticeColorMatrix f3;
    LatticeColorMatrix f4;
    LatticeColorMatrix f5;
    
    START_CODE();
  
    if ( Nd != 4 )
      QDP_error_exit("expecting Nd == 4", Nd);
  
    if ( Ns != 4 )
      QDP_error_exit("expecting Ns == 4", Ns);
  
            
    /* Multiply in the appropriate clover coefficient */
    switch (param.aniso.t_dir)
    {
    case 1:
      f0 = f[0] * param.clovCoeffT;
      f1 = f[1] * param.clovCoeffR;
      f2 = f[2] * param.clovCoeffR;
      f3 = f[3] * param.clovCoeffT;
      f4 = f[4] * param.clovCoeffT;
      f5 = f[5] * param.clovCoeffR;
      break;

    case 2:
      f0 = f[0] * param.clovCoeffR;
      f1 = f[1] * param.clovCoeffT;
      f2 = f[2] * param.clovCoeffR;
      f3 = f[3] * param.clovCoeffT;
      f4 = f[4] * param.clovCoeffR;
      f5 = f[5] * param.clovCoeffT;
      break;

    case 3:
      f0 = f[0] * param.clovCoeffR;
      f1 = f[1] * param.clovCoeffR;
      f2 = f[2] * param.clovCoeffT;
      f3 = f[3] * param.clovCoeffR;
      f4 = f[4] * param.clovCoeffT;
      f5 = f[5] * param.clovCoeffT;
      break;

    default:
      QDP_error_exit("invalid time direction", t_dir);
    }


    /*# Construct diagonal */
    for(int site; site < QDP::Layout::sitesOnNode(); ++site)
    {
      for(int ii = 0; ii < 2*Nc; ii++)
	for(int jj = 0; jj < 2; jj++)
	{
	  tri[site].diag[jj][ii] = 1.0;
	}
    }

    /* The appropriate clover coeffients are already included in the
       field strengths F(mu,nu)! */
    for(int site; site < QDP::Layout::sitesOnNode(); ++site)
    {
      InnerComplex E_minus;
      InnerComplex B_minus;
      InnerComplex ctmp_0;
      InnerComplex ctmp_1;
      InnerReal rtmp_0;
      InnerReal rtmp_1;
      Complex ione;
                              
      FILL(ione,I);
  
      for(int i = 0; i < Nc; ++i)
      {
	/*# diag_L(i,0) = 1 - i*diag(E_z - B_z) */
	/*#             = 1 - i*diag(F(3,2) - F(1,0)) */
	ctmp_0 = _SP_f5[i][i];
	ctmp_0 -= _SP_f0[i][i];
	rtmp_0 = imag(ctmp_0);
	_DIAG_SP_tri[0][i] += rtmp_0;

	/*# diag_L(i+Nc,0) = 1 + i*diag(E_z - B_z) */
	/*#                = 1 + i*diag(F(3,2) - F(1,0)) */
	_DIAG_SP_tri[0][i+Nc] -= rtmp_0;

	/*# diag_L(i,1) = 1 + i*diag(E_z + B_z) */
	/*#             = 1 + i*diag(F(3,2) + F(1,0)) */
	ctmp_1 = _SP_f5[i][i];
	ctmp_1 += _SP_f0[i][i];
	rtmp_1 = imag(ctmp_1);
	_DIAG_SP_tri[1][i] -= rtmp_1;

	/*# diag_L(i+Nc,1) = 1 - i*diag(E_z + B_z) */
	/*#                = 1 - i*diag(F(3,2) + F(1,0)) */
	_DIAG_SP_tri[1][i+Nc] += rtmp_1;
      }

      /*# Construct lower triangular portion */
      /*# Block diagonal terms */
      for(int i = 1; i < Nc; ++i)
      {
	for(int j = 0; j < i; ++j)
	{
	  int elem_ij  = i*(i-1)/2 + j;
	  int elem_tmp = (i+Nc)*(i+Nc-1)/2 + j+Nc;

	  /*# L(i,j,0) = -i*(E_z - B_z)[i,j] */
	  /*#          = -i*(F(3,2) - F(1,0)) */
	  ctmp_0 = _SP_f0[j][i];
	  ctmp_0 -= _SP_f5[j][i];
	  _OFFD_SP_tri[0][elem_ij] = ctmp_0 * ione;

	  /*# L(i+Nc,j+Nc,0) = +i*(E_z - B_z)[i,j] */
	  /*#                = +i*(F(3,2) - F(1,0)) */
	  _OFFD_SP_tri[0][elem_tmp] = -_OFFD_SP_tri[0][elem_ij];

	  /*# L(i,j,1) = i*(E_z + B_z)[i,j] */
	  /*#          = i*(F(3,2) + F(1,0)) */
	  ctmp_1 = _SP_f5[j][i];
	  ctmp_1 += _SP_f0[j][i];
	  _OFFD_SP_tri[1][elem_ij] = ctmp_1 * ione;

	  /*# L(i+Nc,j+Nc,1) = -i*(E_z + B_z)[i,j] */
	  /*#                = -i*(F(3,2) + F(1,0)) */
	  _OFFD_SP_tri[1][elem_tmp] = -_OFFD_SP_tri[1][elem_ij];
	}
      }
  
      /*# Off-diagonal */
      for(int i = 0; i < Nc; ++i)
      {
	for(int j = 0; j < Nc; ++j)
	{
	  int elem_ij  = (i+Nc)*(i+Nc-1)/2 + j;

	  /*# i*E_- = (i*E_x + E_y) */
	  /*#       = (i*F(3,0) + F(3,1)) */
	  E_minus = _SP_f2[j][i] * ione;
	  E_minus += _SP_f4[j][i];

	  /*# i*B_- = (i*B_x + B_y) */
	  /*#       = (i*F(2,1) - F(2,0)) */
	  B_minus = _SP_f3[j][i] * ione;
	  B_minus -= _SP_f1[j][i];

	  /*# L(i+Nc,j,0) = -i*(E_- - B_-)  */
	  _OFFD_SP_tri[0][elem_ij] = B_minus - E_minus;

	  /*# L(i+Nc,j,1) = +i*(E_- + B_-)  */
	  _OFFD_SP_tri[1][elem_ij] = E_minus + B_minus;
	}
      }
    }
              
    END_CODE();
  }
  

  //! Invert
  /*!
   * Computes the inverse of the term on cb using Cholesky
   */
  void choles(int cb)
  {
    Double logdet;
    chlclovms(false, logdet, cb);
  }


  //! Invert
  /*!
   * Computes the inverse of the term on cb using Cholesky
   *
   * \return logarithm of the determinant  
   */
  Double cholesDet(int cb)
  {
    Double logdet;
    chlclovms(true, logdet, cb);
    return det;
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
  void QDPCloverTerm::choles(bool DetP, Double& logdet, int cb)
  {
    LatticeReal log_diag;

    START_CODE();

    if ( 2*Nc < 3 )
      QDP_error_exit("Matrix is too small", Nc, Ns);
  
    log_diag = 0;
  
    int n = 2*Nc;

    /*# Cholesky decompose  A = L.L^dag */
    /*# NOTE!!: I can store this matrix in  invclov, but will need a */
    /*#   temporary  diag */
    for(int site=0; site < QDP::Layout::sitesOnNode(); ++site)
    {
      multi1d<InnerReal> diag_g(n);
      multi1d<InnerComplex> v1(n);
      InnerComplex sum;
      InnerReal one;
      InnerReal zero;
      InnerReal lrtmp;

      one = 1;
      zero = 0;
  
      for(int s = 0; s < 2; ++s)
      {
	int elem_jk = 0;
    
	for(int j = 0; j <  n; ++j)
	{
	  /*# Multiply clover mass term against basis vector.  */
	  /*# Actually, I need a column of the lower triang matrix clov. */
	  v1[j] = cmplx(_DIAG_SP_clo[s][j],zero);
    
	  int elem_ij = elem_jk + 2*j;
      
	  for(int i = j+1; i < n; ++i)
	  {
	    v1[i] = _OFFD_SP_clo[s][elem_ij];
	    elem_ij += i;
	  }
      
	  /*# Back to cholesky */
	  /*# forward substitute */
	  for(int k = 0; k < j; ++k)
	  {
	    int elem_ik = elem_jk;
	
	    for(int i = j; i < n; ++i)
	    {
	      v1[i] -= adj(_OFFD_SP_invcl[s][elem_jk]) * _OFFD_SP_invcl[s][elem_ik];
	      elem_ik += i;
	    }
	    elem_jk++;
	  }

	  /*# The diagonal is (should be!!) real and positive */
	  diag_g[j] = real(v1[j]);
	  
	  /*#+ */
	  /*# Squeeze in computation of log(Det) */
	  /*#- */
	  lrtmp = log(diag_g[j]);
	  _SP_log_diag += lrtmp;
	  
	  diag_g[j] = sqrt(diag_g[j]);
	  diag_g[j] = one / diag_g[j];
      
	  /*# backward substitute */
	  int elem_ij = elem_jk + j;
	  for(int i = j+1; i < n; ++i)
	  {
	    _OFFD_SP_invcl[s][elem_ij] = v1[i] * diag_g[j];
	    elem_ij += i;
	  }
	}

	/*# Use forward and back substitution to construct  _OFFD_SP_invcl = lower(A^-1) */
	for(int k = 0; k < n; ++k)
	{
	  for(int i = 0; i < k; ++i)
	    v1[i] = 0;
	  
	  /*# Forward substitution */
	  v1[k] = cmplx(diag_g[k],zero);
      
	  for(int i = k+1; i < n; ++i)
	  {
	    sum = 0;
	    elem_ij = i*(i-1)/2+k;	
	    for(j = k; j < i; ++j)
	    {
	      sum -= _OFFD_SP_invcl[s][elem_ij] * v1[j];
	      elem_ij++;
	    }
	
	    v1[i] = sum * diag_g[i];
	  }
      
	  /*# Backward substitution */
	  v1[n-1] = v1[n-1] * diag_g[n-1];
     
	  for(int i = n-2; (int)i >= (int)k; --i)
	  {
	    sum = v1[i];
	    
	    int elem_ji = ((i+1)*i)/2+i;
	    for(int j = i+1; j < n; ++j)
	    {
	      sum -= adj(_OFFD_SP_invcl[s][elem_ji]) * v1[j];
	      elem_ji += j;
	    }
	    v1[i] = sum * diag_g[i];
	  }

	  /*# Overwrite column k of _OFFD_SP_invcl */
	  _DIAG_SP_invcl[s][k] = real(v1[k]);

	  elem_ik = ((k+1)*k)/2+k;
      
	  for(i = k+1; i < n; ++i)
	  {
	    _OFFD_SP_invcl[s][elem_ik] = v1[i];
	    elem_ik += i;
	  }
	}
      }
    }
  
    if ( DetP )
    {
      logdet = sum(log_diag);
    }
    else
    {
      logdet = 0;
    }
  
    END_CODE();
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
   * \param cb      Checkerboard of OUTPUT vector               (Read) 
   */
  void QDPCloverTerm::apply(LatticeFermion& chi, const LatticeFermion& psi, 
			    enum PlusMinus isign, int cb) const
  {
    START_CODE();

    if ( Ns != 4 )
      QDP_error_exit("code requires Ns == 4", Ns);

    int n = 2*Nc;

    for(site; site < QDP::Layout::sitesOnNode(); ++site)
    {
      for(int i = 0; i < n; ++i)
      {
	static_cast<PrimitiveLDUFerm&>(chi.elem(site))[0][i] = _DIAG_SP_clov[0][i] * _SP_psi[0][i];
	chi.elem(site).elem[2][i] = _DIAG_SP_clov[1][i] * _SP_psi[2][i];
      }

      int kij=0;  
      for(int i = 0; i < n; ++i)
      {
	for(int j=0; j < i; j++)
	{
	  _SP_chi[0][i] += _OFFD_SP_clov[0][kij] * _SP_psi[0][j];
	  _SP_chi[0][j] += adj(_OFFD_SP_clov[0][kij]) * _SP_psi[0][i];
	  _SP_chi[2][i] += _OFFD_SP_clov[1][kij] * _SP_psi[2][j];
	  _SP_chi[2][j] += adj(_OFFD_SP_clov[1][kij]) * _SP_psi[2][i];
	  kij++;
	}
      }
    }


    END_CODE();
  }

}
