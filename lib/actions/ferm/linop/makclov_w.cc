#error "NOT FULLY CONVERTED"

/*
 * $Id: makclov_w.cc,v 1.1 2005-02-17 02:52:36 edwards Exp $

 * MAKCLOV - calculates

 *    1 - (1/4)*Kappa*sigma(mu,nu) F(mu,nu)

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

 *  sigF(0) = i*Kappa*(E_z + B_z)
 *          = i*Kappa*(F(2,3) + F(0,1))
 *  sigF(1) = i*Kappa*(E_z - B_z)
 *          = i*Kappa*(F(2,3) - F(0,1))
 *  sigF(2) = i*Kappa*(E_+ + B_+)
 *  sigF(3) = i*Kappa*(E_+ - B_+)
 *  sigF(4) = i*Kappa*(E_- + B_-)
 *  sigF(5) = i*Kappa*(E_- - B_-)
 *  i*Kappa*E_+ = Kappa*(i*E_x - E_y)
 *              = Kappa*(i*F(0,3) - F(1,3))
 *  i*Kappa*E_- = Kappa*(i*E_x + E_y)
 *              = Kappa*(i*F(0,3) + F(1,3))
 *  i*Kappa*B_+ = Kappa*(i*B_x - B_y)
 *              = Kappa*(i*F(1,2) + F(0,2))
 *  i*Kappa*B_- = Kappa*(i*B_x + B_y)
 *              = Kappa*(i*F(1,2) - F(0,2))

 *  NOTE: I am using  i*F  of the usual F defined by UKQCD, Heatlie et.al.

 * Arguments:

 *  f      -- field strength tensor F(cb,mu,nu) (Read)
 *  L      -- lower((1/4)*F(cb,mu,nu) sigma(mu,nu))    (Write)
 *  diag_L -- diag(L)                           (Write)
 *  Kappa  -- yep, its kappa                    (Read)
 *  cb     -- checkerboard                      (Read)
 */

include(types.mh)

SUBROUTINE(makclov, f, L, diag_L, Kappa, cb)

LATTICE_FIELD_STRENGTH(f);
LATTICE_TRIANGULAR(L);
LATTICE_DIAG_TRIANGULAR(diag_L);
Real Kappa;
int cb;
{				/* Local variables */
  include(COMMON_DECLARATIONS)

  LatticeColorMatrix fx;
  LatticeSpinMatrix f0;
  LatticeSpinMatrix f1;
  LatticeSpinMatrix f2;
  LatticeSpinMatrix f3;
  LatticeSpinMatrix f4;
  LatticeSpinMatrix f5;
  LatticeComplex E_minus;
  LatticeComplex B_minus;
  LatticeComplex ctmp_0;
  LatticeReal rtmp_0;
  Complex iKappa;
  Real zero;
  int i;
  int j;
  int elem_ij;
  int elem_ji;
  int elem_tmp;
  int n;

  START_CODE();
  
  if ( Nd != 4 )
    QDP_error_exit("expecting Nd == 4", Nd);
  
  if ( Ns != 4 )
    QDP_error_exit("expecting Ns == 4", Ns);
  
  n = 2*Nc;
  
  diag_L = 1;
  
              
  fx = f[0][cb];
  f0 = CAST(fx);
  fx = f[1][cb];
  f1 = CAST(fx);
  fx = f[2][cb];
  f2 = CAST(fx);
  fx = f[3][cb];
  f3 = CAST(fx);
  fx = f[4][cb];
  f4 = CAST(fx);
  fx = f[5][cb];
  f5 = CAST(fx);

              
  zero = 0;
  iKappa = cmplx(zero,Kappa);
  
  /*# Construct diagonal */
  for(i = 0; i < Nc; ++i)
  {
    /*# diag_L(i,0) = 1 - i*Kappa*diag(E_z - B_z) */
    /*#             = 1 - i*Kappa*diag(F(2,3) - F(0,1)) */
    ctmp_0 = f5[i][i] - f0[i][i];
    rtmp_0 = imag(ctmp_0);
    diag_L[0][i] += rtmp_0 * Kappa;

    /*# diag_L(i+Nc,0) = 1 + i*Kappa*diag(E_z - B_z) */
    /*#                = 1 + i*Kappa*diag(F(2,3) - F(0,1)) */
    diag_L[0][i+Nc] -= rtmp_0 * Kappa;

    /*# diag_L(i,1) = 1 + i*Kappa*diag(E_z + B_z) */
    /*#             = 1 + i*Kappa*diag(F(2,3) + F(0,1)) */
    ctmp_0 = f5[i][i] + f0[i][i];
    rtmp_0 = imag(ctmp_0);
    diag_L[1][i] -= rtmp_0 * Kappa;

    /*# diag_L(i+Nc,1) = 1 - i*Kappa*diag(E_z + B_z) */
    /*#                = 1 - i*Kappa*diag(F(2,3) + F(0,1)) */
    diag_L[1][i+Nc] += rtmp_0 * Kappa;
  }
  
  /*# Construct lower triangular portion */
  /*# Block diagonal terms */
  for(i = 1; i < Nc; ++i)
  {
    for(j = 0; j < i; ++j)
    {
      elem_ij  = i*(i-1)/2 + j;
      elem_tmp = (i+Nc)*(i+Nc-1)/2 + j+Nc;

      /*# L(i,j,0) = -i*Kappa*(E_z - B_z)[i,j] */
      /*#          = -i*Kappa*(F(2,3) - F(0,1)) */
      ctmp_0 = f0[j][i] - f5[j][i];
      L[0][elem_ij] = ctmp_0 * iKappa;

      /*# L(i+Nc,j+Nc,0) = +i*Kappa*(E_z - B_z)[i,j] */
      /*#                = +i*Kappa*(F(2,3) - F(0,1)) */
      L[0][elem_tmp] = -L[0][elem_ij];

      /*# L(i,j,1) = i*Kappa*(E_z + B_z)[i,j] */
      /*#          = i*Kappa*(F(2,3) + F(0,1)) */
      ctmp_0 = f0[j][i] + f5[j][i];
      L[1][elem_ij] = ctmp_0 * iKappa;

      /*# L(i+Nc,j+Nc,1) = -i*Kappa*(E_z + B_z)[i,j] */
      /*#                = -i*Kappa*(F(2,3) + F(0,1)) */
      L[1][elem_tmp] = -L[1][elem_ij];
    }
  }
  
  /*# Off-diagonal */
  for(i = 0; i < Nc; ++i)
  {
    for(j = 0; j < Nc; ++j)
    {
      elem_ij  = (i+Nc)*(i+Nc-1)/2 + j;

      /*# i*Kappa*E_- = Kappa*(i*E_x + E_y) */
      /*#             = Kappa*(i*F(0,3) + F(1,3)) */
      E_minus = f2[j][i] * iKappa;
      E_minus += f4[j][i] * Kappa;

      /*# i*Kappa*B_- = Kappa*(i*B_x + B_y) */
      /*#             = Kappa*(i*F(1,2) - F(0,2)) */
      B_minus = f3[j][i] * iKappa;
      B_minus -= f1[j][i] * Kappa;

      /*# L(i+Nc,j,0) = -i*Kappa*(E_- - B_-)  */
      L[0][elem_ij] = B_minus - E_minus;

      /*# L(i+Nc,j,1) = +i*Kappa*(E_- + B_-)  */
      L[1][elem_ij] = E_minus + B_minus;
    }
  }
  
                        
  END_CODE();
}
