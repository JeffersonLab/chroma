// $Id: clover_term_qdp_w.cc,v 2.4 2006-01-01 05:12:30 edwards Exp $
/*! \file
 *  \brief Clover term linear operator
 *
 *  This particular implementation is specific to a scalar-like
 *  architecture
 */

#include "chromabase.h"
#include "actions/ferm/linop/clover_term_qdp_w.h"
#include "meas/glue/mesfield.h"


namespace Chroma 
{ 


  //! XML output
  inline
  XMLWriter& operator<<(XMLWriter& xml, const PrimitiveClovTriang& d)
  {
    xml.openTag("PrimClovTriang");

    XMLWriterAPI::AttributeList alist;

    xml.openTag("Diag");
    for(int i=0; i < 2; ++i)
    {
      for(int j=0; j < 2*Nc; ++j)
      {
	alist.clear();
	alist.push_back(XMLWriterAPI::Attribute("block", i));
	alist.push_back(XMLWriterAPI::Attribute("col", j));

	xml.openTag("elem", alist);
	xml << d.diag[i][j];
	xml.closeTag();
      }
    }
    xml.closeTag(); // Diag

    xml.openTag("Offd");
    for(int i=0; i < 2; ++i)
    {
      for(int j=0; j < 2*Nc*Nc-Nc; ++j)
      {
	alist.clear();
	alist.push_back(XMLWriterAPI::Attribute("block", i));
	alist.push_back(XMLWriterAPI::Attribute("col", j));

	xml.openTag("elem", alist);
	xml << d.offd[i][j];
	xml.closeTag();
      }
    }
    xml.closeTag(); // Offd

    xml.closeTag(); // PrimClovTriang
    return xml;
  }


  // Reader/writers
  void read(XMLReader& xml, const string& path, PrimitiveClovTriang& param)
  {
    QDP_error_exit("clover reader not implemented");
  }

  void write(XMLWriter& xml, const string& path, const PrimitiveClovTriang& param)
  {
    push(xml,path);
    xml << param;
    pop(xml);
  }



  //! Creation routine
  void QDPCloverTerm::create(const multi1d<LatticeColorMatrix>& u_, 	
			     const CloverFermActParams& param_)
  {
    u = u_;
    param = param_;

    //
    // Yuk. Some bits of knowledge of the dslash term are buried in the 
    // effective mass term. They show up here. If I wanted some more 
    // complicated dslash then this will have to be fixed/adjusted.
    //
    Real diag_mass;
    {
      Real ff = where(param.anisoParam.anisoP, param.anisoParam.nu / param.anisoParam.xi_0, Real(1));
      diag_mass = 1 + (Nd-1)*ff + param.Mass;
    }

    /* Calculate F(mu,nu) */
    multi1d<LatticeColorMatrix> f;
    mesField(f, u_);
    makeClov(f, diag_mass);
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
   *  \param diag_mass effective mass term                      (Read)
   */
  void QDPCloverTerm::makeClov(const multi1d<LatticeColorMatrix>& f, const Real& diag_mass)
  {
//    QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << endl;

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
    switch (param.anisoParam.t_dir)
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
      QDPIO::cerr << __func__ << ": invalid time direction: t_dir= " 
		  << param.anisoParam.t_dir << endl;
      QDP_abort(1);
    }


    tri.resize(QDP::Layout::sitesOnNode());  // hold local lattice

    /*# Construct diagonal */
    for(int site = 0; site < QDP::Layout::sitesOnNode(); ++site)
    {
      for(int jj = 0; jj < 2; jj++)
      {
	for(int ii = 0; ii < 2*Nc; ii++)
	{
	  tri[site].diag[jj][ii] = diag_mass.elem().elem().elem();
	}
      }
    }


    /* The appropriate clover coeffients are already included in the
       field strengths F(mu,nu)! */
    for(int site = 0; site < QDP::Layout::sitesOnNode(); ++site)
    {
      RComplex<REAL> E_minus;
      RComplex<REAL> B_minus;
      RComplex<REAL> ctmp_0;
      RComplex<REAL> ctmp_1;
      RScalar<REAL> rtmp_0;
      RScalar<REAL> rtmp_1;
  
      for(int i = 0; i < Nc; ++i)
      {
	/*# diag_L(i,0) = 1 - i*diag(E_z - B_z) */
	/*#             = 1 - i*diag(F(3,2) - F(1,0)) */
	ctmp_0 = f5.elem(site).elem().elem(i,i);
	ctmp_0 -= f0.elem(site).elem().elem(i,i);
	rtmp_0 = imag(ctmp_0);
	tri[site].diag[0][i] += rtmp_0;

	/*# diag_L(i+Nc,0) = 1 + i*diag(E_z - B_z) */
	/*#                = 1 + i*diag(F(3,2) - F(1,0)) */
	tri[site].diag[0][i+Nc] -= rtmp_0;

	/*# diag_L(i,1) = 1 + i*diag(E_z + B_z) */
	/*#             = 1 + i*diag(F(3,2) + F(1,0)) */
	ctmp_1 = f5.elem(site).elem().elem(i,i);
	ctmp_1 += f0.elem(site).elem().elem(i,i);
	rtmp_1 = imag(ctmp_1);
	tri[site].diag[1][i] -= rtmp_1;

	/*# diag_L(i+Nc,1) = 1 - i*diag(E_z + B_z) */
	/*#                = 1 - i*diag(F(3,2) + F(1,0)) */
	tri[site].diag[1][i+Nc] += rtmp_1;
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
	  ctmp_0 = f0.elem(site).elem().elem(j,i);
	  ctmp_0 -= f5.elem(site).elem().elem(j,i);
	  tri[site].offd[0][elem_ij] = timesI(ctmp_0);

	  /*# L(i+Nc,j+Nc,0) = +i*(E_z - B_z)[i,j] */
	  /*#                = +i*(F(3,2) - F(1,0)) */
	  tri[site].offd[0][elem_tmp] = -tri[site].offd[0][elem_ij];

	  /*# L(i,j,1) = i*(E_z + B_z)[i,j] */
	  /*#          = i*(F(3,2) + F(1,0)) */
	  ctmp_1 = f5.elem(site).elem().elem(j,i);
	  ctmp_1 += f0.elem(site).elem().elem(j,i);
	  tri[site].offd[1][elem_ij] = timesI(ctmp_1);

	  /*# L(i+Nc,j+Nc,1) = -i*(E_z + B_z)[i,j] */
	  /*#                = -i*(F(3,2) + F(1,0)) */
	  tri[site].offd[1][elem_tmp] = -tri[site].offd[1][elem_ij];
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
	  E_minus = timesI(f2.elem(site).elem().elem(j,i));
	  E_minus += f4.elem(site).elem().elem(j,i);

	  /*# i*B_- = (i*B_x + B_y) */
	  /*#       = (i*F(2,1) - F(2,0)) */
	  B_minus = timesI(f3.elem(site).elem().elem(j,i));
	  B_minus -= f1.elem(site).elem().elem(j,i);

	  /*# L(i+Nc,j,0) = -i*(E_- - B_-)  */
	  tri[site].offd[0][elem_ij] = B_minus - E_minus;

	  /*# L(i+Nc,j,1) = +i*(E_- + B_-)  */
	  tri[site].offd[1][elem_ij] = E_minus + B_minus;
	}
      }
    }
              
    END_CODE();

#if 0
    {
      XMLFileWriter xml("f.xml");
      push(xml,"f");
      write(xml, "clovT", param.clovCoeffT);
      write(xml, "clovR", param.clovCoeffR);
      write(xml,"f",f);
      pop(xml);
    }
    {
      XMLFileWriter xml("makeclov.xml");
      push(xml,"makeclov");
      write(xml,"tri",tri);
      pop(xml);
    }
#endif

//    QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << endl;
  }
  

  //! Invert
  /*!
   * Computes the inverse of the term on cb using Cholesky
   */
  void QDPCloverTerm::choles(int cb)
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
  Double QDPCloverTerm::cholesDet(int cb)
  {
    Double logdet;
    chlclovms(true, logdet, cb);
    return logdet;
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
  void QDPCloverTerm::chlclovms(bool DetP, Double& logdet, int cb)
  {
//    QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << endl;

    START_CODE();

    LatticeReal log_diag;

    if ( 2*Nc < 3 )
      QDP_error_exit("Matrix is too small", Nc, Ns);
  
    log_diag = zero;
  
    int n = 2*Nc;

    /*# Cholesky decompose  A = L.L^dag */
    /*# NOTE!!: I can store this matrix in  invclov, but will need a */
    /*#   temporary  diag */
    for(int ssite=0; ssite < rb[cb].numSiteTable(); ++ssite) 
    {
      int site = rb[cb].siteTable()[ssite];

      PrimitiveClovTriang  invcl;

      multi1d< RScalar<REAL> > diag_g(n);
      multi1d< RComplex<REAL> > v1(n);
      RComplex<REAL> sum;
      RScalar<REAL> one;
      RScalar<REAL> zero;
      RScalar<REAL> lrtmp;

      one = 1;
      zero = 0;
  
      for(int s = 0; s < 2; ++s)
      {
	int elem_jk = 0;
	int elem_ij;

	for(int j = 0; j <  n; ++j)
	{
	  /*# Multiply clover mass term against basis vector.  */
	  /*# Actually, I need a column of the lower triang matrix clov. */
	  v1[j] = cmplx(tri[site].diag[s][j],zero);
    
	  elem_ij = elem_jk + 2*j;
      
	  for(int i = j+1; i < n; ++i)
	  {
	    v1[i] = tri[site].offd[s][elem_ij];
	    elem_ij += i;
	  }
      
	  /*# Back to cholesky */
	  /*# forward substitute */
	  for(int k = 0; k < j; ++k)
	  {
	    int elem_ik = elem_jk;
	
	    for(int i = j; i < n; ++i)
	    {
	      v1[i] -= adj(invcl.offd[s][elem_jk]) * invcl.offd[s][elem_ik];
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
	  log_diag.elem(site).elem().elem() += lrtmp;
	  
	  diag_g[j] = sqrt(diag_g[j]);
	  diag_g[j] = one / diag_g[j];
      
	  /*# backward substitute */
	  elem_ij = elem_jk + j;
	  for(int i = j+1; i < n; ++i)
	  {
	    invcl.offd[s][elem_ij] = v1[i] * diag_g[j];
	    elem_ij += i;
	  }
	}

	/*# Use forward and back substitution to construct  invcl.offd = lower(A^-1) */
	for(int k = 0; k < n; ++k)
	{
	  for(int i = 0; i < k; ++i)
	    zero_rep(v1[i]);
	  
	  /*# Forward substitution */
	  v1[k] = cmplx(diag_g[k],zero);
      
	  for(int i = k+1; i < n; ++i)
	  {
	    zero_rep(sum);
	    elem_ij = i*(i-1)/2+k;	
	    for(int j = k; j < i; ++j)
	    {
	      sum -= invcl.offd[s][elem_ij] * v1[j];
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
	      sum -= adj(invcl.offd[s][elem_ji]) * v1[j];
	      elem_ji += j;
	    }
	    v1[i] = sum * diag_g[i];
	  }

	  /*# Overwrite column k of invcl.offd */
	  invcl.diag[s][k] = real(v1[k]);

	  int elem_ik = ((k+1)*k)/2+k;
      
	  for(int i = k+1; i < n; ++i)
	  {
	    invcl.offd[s][elem_ik] = v1[i];
	    elem_ik += i;
	  }
	}
      }

      // Overwrite original element
      tri[site] = invcl;
    }
  
    // Overwrite

    if ( DetP )
    {
      logdet = sum(log_diag);
    }
    else
    {
      logdet = 0;
    }
  
    END_CODE();

//    QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << endl;
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

    const multi1d<int>& tab = rb[cb].siteTable();
    for(int ssite=0; ssite < tab.size(); ++ssite) 
    {
      int site = tab[ssite];

      RComplex<REAL>* cchi = (RComplex<REAL>*)&(chi.elem(site).elem(0).elem(0));
      const RComplex<REAL>* ppsi = (const RComplex<REAL>*)&(psi.elem(site).elem(0).elem(0));

      for(int i = 0; i < n; ++i)
      {
	cchi[0*n+i] = tri[site].diag[0][i] * ppsi[0*n+i];
	cchi[1*n+i] = tri[site].diag[1][i] * ppsi[1*n+i];
      }

      int kij = 0;  
      for(int i = 0; i < n; ++i)
      {
	for(int j = 0; j < i; j++)
	{
	  cchi[0*n+i] += tri[site].offd[0][kij] * ppsi[0*n+j];
	  cchi[0*n+j] += adj(tri[site].offd[0][kij]) * ppsi[0*n+i];
	  cchi[1*n+i] += tri[site].offd[1][kij] * ppsi[1*n+j];
	  cchi[1*n+j] += adj(tri[site].offd[1][kij]) * ppsi[1*n+i];
	  kij++;
	}
      }
    }


    END_CODE();
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
  void QDPCloverTerm::triacntr(LatticeColorMatrix& B, int mat, int cb) const
  {
    START_CODE();
  
    if ( mat < 0  ||  mat > 15 )
    {
      QDPIO::cerr << __func__ << ": Gamma out of range: mat = " << mat << endl;
      QDP_abort(1);
    }
  
    switch( mat )
    {
    case 0:
      /*# gamma(   0)   1  0  0  0            # ( 0000 )  --> 0 */
      /*#               0  1  0  0 */
      /*#               0  0  1  0 */
      /*#               0  0  0  1 */
      /*# From diagonal part */
      for(int ssite=0; ssite < rb[cb].numSiteTable(); ++ssite) 
      {
	int site = rb[cb].siteTable()[ssite];

	RComplex<REAL> lctmp0;
	RScalar<REAL> lr_zero0;
	RScalar<REAL> lrtmp0;
  
	lr_zero0 = 0;
  
	for(int i0 = 0; i0 < Nc; ++i0)
	{
	  lrtmp0 = tri[site].diag[0][i0];
	  lrtmp0 += tri[site].diag[0][i0+Nc];
	  lrtmp0 += tri[site].diag[1][i0];
	  lrtmp0 += tri[site].diag[1][i0+Nc];
	  B.elem(site).elem().elem(i0,i0) = cmplx(lrtmp0,lr_zero0);
	}

	/*# From lower triangular portion */
	int elem_ij0 = 0;
	for(int i0 = 1; i0 < Nc; ++i0)
	{
	  int elem_ijb0 = (i0+Nc)*(i0+Nc-1)/2 + Nc;

	  for(int j0 = 0; j0 < i0; ++j0)
	  {
	    lctmp0 = tri[site].offd[0][elem_ij0];
	    lctmp0 += tri[site].offd[0][elem_ijb0];
	    lctmp0 += tri[site].offd[1][elem_ij0];
	    lctmp0 += tri[site].offd[1][elem_ijb0];

	    B.elem(site).elem().elem(j0,i0) = lctmp0;
	    B.elem(site).elem().elem(i0,j0) = adj(lctmp0);
	      
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
      for(int ssite=0; ssite < rb[cb].numSiteTable(); ++ssite) 
      {
	int site = rb[cb].siteTable()[ssite];

	RComplex<REAL> lctmp3;
	RScalar<REAL> lr_zero3;
	RScalar<REAL> lrtmp3;
                          
	lr_zero3 = 0;
  
	for(int i3 = 0; i3 < Nc; ++i3)
	{
	  lrtmp3 = tri[site].diag[0][i3+Nc];
	  lrtmp3 -= tri[site].diag[0][i3];
	  lrtmp3 -= tri[site].diag[1][i3];
	  lrtmp3 += tri[site].diag[1][i3+Nc];
	  B.elem(site).elem().elem(i3,i3) = cmplx(lr_zero3,lrtmp3);
	}
	
	/*# From lower triangular portion */
	int elem_ij3 = 0;
	for(int i3 = 1; i3 < Nc; ++i3)
	{
	  int elem_ijb3 = (i3+Nc)*(i3+Nc-1)/2 + Nc;

	  for(int j3 = 0; j3 < i3; ++j3)
	  {
	    lctmp3 = tri[site].offd[0][elem_ijb3];
	    lctmp3 -= tri[site].offd[0][elem_ij3];
	    lctmp3 -= tri[site].offd[1][elem_ij3];
	    lctmp3 += tri[site].offd[1][elem_ijb3];

	    B.elem(site).elem().elem(j3,i3) = timesI(lctmp3);
	    B.elem(site).elem().elem(i3,j3) = timesI(adj(lctmp3));
	    
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
      for(int ssite=0; ssite < rb[cb].numSiteTable(); ++ssite) 
      {
	int site = rb[cb].siteTable()[ssite];

	RComplex<REAL> lctmp5;
	RScalar<REAL> lrtmp5;
                          
	for(int i5 = 0; i5 < Nc; ++i5)
	{
	  int elem_ij5 = (i5+Nc)*(i5+Nc-1)/2;

	  for(int j5 = 0; j5 < Nc; ++j5)
	  {
	    int elem_ji5 = (j5+Nc)*(j5+Nc-1)/2 + i5;
	  
	    lctmp5 = adj(tri[site].offd[0][elem_ji5]);
	    lctmp5 -= tri[site].offd[0][elem_ij5];
	    lctmp5 += adj(tri[site].offd[1][elem_ji5]);
	    lctmp5 -= tri[site].offd[1][elem_ij5];

	    B.elem(site).elem().elem(j5,i5) = lctmp5;

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
      for(int ssite=0; ssite < rb[cb].numSiteTable(); ++ssite) 
      {
	int site = rb[cb].siteTable()[ssite];

	RComplex<REAL> lctmp6;
	RScalar<REAL> lrtmp6;
                          
	for(int i6 = 0; i6 < Nc; ++i6)
	{
	  int elem_ij6 = (i6+Nc)*(i6+Nc-1)/2;

	  for(int j6 = 0; j6 < Nc; ++j6)
	  {
	    int elem_ji6 = (j6+Nc)*(j6+Nc-1)/2 + i6;

	    lctmp6 = adj(tri[site].offd[0][elem_ji6]);
	    lctmp6 += tri[site].offd[0][elem_ij6];
	    lctmp6 += adj(tri[site].offd[1][elem_ji6]);
	    lctmp6 += tri[site].offd[1][elem_ij6];

	    B.elem(site).elem().elem(j6,i6) = timesMinusI(lctmp6);

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
      for(int ssite=0; ssite < rb[cb].numSiteTable(); ++ssite) 
      {
	int site = rb[cb].siteTable()[ssite];

	RComplex<REAL> lctmp9;
	RScalar<REAL> lrtmp9;
                          
	for(int i9 = 0; i9 < Nc; ++i9)
	{
	  int elem_ij9 = (i9+Nc)*(i9+Nc-1)/2;

	  for(int j9 = 0; j9 < Nc; ++j9)
	  {
	    int elem_ji9 = (j9+Nc)*(j9+Nc-1)/2 + i9;

	    lctmp9 = adj(tri[site].offd[0][elem_ji9]);
	    lctmp9 += tri[site].offd[0][elem_ij9];
	    lctmp9 -= adj(tri[site].offd[1][elem_ji9]);
	    lctmp9 -= tri[site].offd[1][elem_ij9];

	    B.elem(site).elem().elem(j9,i9) = timesI(lctmp9);

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
      for(int ssite=0; ssite < rb[cb].numSiteTable(); ++ssite) 
      {
	int site = rb[cb].siteTable()[ssite];

	RComplex<REAL> lctmp10;
	RScalar<REAL> lrtmp10;
                          
	for(int i10 = 0; i10 < Nc; ++i10)
	{
	  int elem_ij10 = (i10+Nc)*(i10+Nc-1)/2;
	
	  for(int j10 = 0; j10 < Nc; ++j10)
	  {
	    int elem_ji10 = (j10+Nc)*(j10+Nc-1)/2 + i10;

	    lctmp10 = adj(tri[site].offd[0][elem_ji10]);
	    lctmp10 -= tri[site].offd[0][elem_ij10];
	    lctmp10 -= adj(tri[site].offd[1][elem_ji10]);
	    lctmp10 += tri[site].offd[1][elem_ij10];

	    B.elem(site).elem().elem(j10,i10) = lctmp10;

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
      for(int ssite=0; ssite < rb[cb].numSiteTable(); ++ssite) 
      {
	int site = rb[cb].siteTable()[ssite];

	RComplex<REAL> lctmp12;
	RScalar<REAL> lr_zero12;
	RScalar<REAL> lrtmp12;
                          
	lr_zero12 = 0;
  
	for(int i12 = 0; i12 < Nc; ++i12)
	{
	  lrtmp12 = tri[site].diag[0][i12];
	  lrtmp12 -= tri[site].diag[0][i12+Nc];
	  lrtmp12 -= tri[site].diag[1][i12];
	  lrtmp12 += tri[site].diag[1][i12+Nc];
	  B.elem(site).elem().elem(i12,i12) = cmplx(lr_zero12,lrtmp12);
	}
    
	/*# From lower triangular portion */
	int elem_ij12 = 0;
	for(int i12 = 1; i12 < Nc; ++i12)
	{
	  int elem_ijb12 = (i12+Nc)*(i12+Nc-1)/2 + Nc;

	  for(int j12 = 0; j12 < i12; ++j12)
	  {
	    lctmp12 = tri[site].offd[0][elem_ij12];
	    lctmp12 -= tri[site].offd[0][elem_ijb12];
	    lctmp12 -= tri[site].offd[1][elem_ij12];
	    lctmp12 += tri[site].offd[1][elem_ijb12];
	
	    B.elem(site).elem().elem(j12,i12) = timesI(lctmp12);
	    B.elem(site).elem().elem(i12,j12) = timesI(adj(lctmp12));
	
	    elem_ij12++;
	    elem_ijb12++;
	  }
	}
      }
      break;
    
    default:
      B = zero;
    }
  

    END_CODE();
  }
}

