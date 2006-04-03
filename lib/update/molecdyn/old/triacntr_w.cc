/*# + */
/*# $Id: triacntr_w.cc,v 3.0 2006-04-03 04:59:10 edwards Exp $ ($Date: 2006-04-03 04:59:10 $) */

/*# This routine is specific to Wilson fermions! */
/*# TRIACNTR - calculates */

/*#     Tr_D ( Gamma_mat L ) */

/*#  the trace over the Dirac indices for one of the 16 Gamma matrices */
/*#  and a hermitian color x spin matrix A, stored as a block diagonal */
/*#  complex lower triangular matrix L and a real diagonal diag_L. */

/*#  Here 0 <= mat <= 15 and */
/*#  if mat = mat_1 + mat_2 * 2 + mat_3 * 4 + mat_4 * 8 */

/*#  Gamma(mat) = gamma(1)^(mat_1) * gamma(2)^(mat_2) * gamma(3)^(mat_3) */
/*#             * gamma(4)^(mat_4) */

/*#  Further, in basis for the Gamma matrices used, A is of the form */

/*#      | A_0 |  0  | */
/*#  A = | --------- | */
/*#      |  0  | A_1 | */


/*# Arguments: */

/*#  L       -- lower triangular matrix		(Read) */
/*#  diag_L  -- diag(L)				(Read) */
/*#  mat     -- label of the Gamma matrix		(Read) */
/*#  B       -- the resulting SU(N) color matrix	(Write) */

void triacntr(LatticeColorMatrix& B,
	      const LATTICE_TRIANG& clov,
	      int mat)
{
  START_CODE("subroutine");;
  
  if ( mat < 0  ||  mat > 15 )
    QDP_error_exit("Gamma out of range", mat);
  
  switch( mat )
  {
  case 0:
    /*# gamma(   0)   1  0  0  0            # ( 0000 )  --> 0 */
    /*#               0  1  0  0 */
    /*#               0  0  1  0 */
    /*#               0  0  0  1 */
    /*# From diagonal part */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp0;
	InnerReal lr_zero0;
	InnerReal lrtmp0;
	Complex sqm0;
	int i0;
	int j0;
	int elem_ij0;
	int elem_ijb0;
  
                          
	lr_zero0 = 0;
	FILL(sqm0, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i0 = 0; i0 < Nc; ++i0)
	    {
	      lrtmp0 = _DIAG_SP_clov[0][i0];
	      lrtmp0 += _DIAG_SP_clov[0][i0+Nc];
	      lrtmp0 += _DIAG_SP_clov[1][i0];
	      lrtmp0 += _DIAG_SP_clov[1][i0+Nc];
	      _SP_B[i0][i0] = cmplx(lrtmp0,lr_zero0);
	    }

	    /*# From lower triangular portion */
	    elem_ij0 = 0;
	    for(i0 = 1; i0 < Nc; ++i0)
	    {
	      elem_ijb0 = (i0+Nc)*(i0+Nc-1)/2 + Nc;
	      for(j0 = 0; j0 < i0; ++j0)
	      {
		lctmp0 = _OFFD_SP_clov[0][elem_ij0];
		lctmp0 += _OFFD_SP_clov[0][elem_ijb0];
		lctmp0 += _OFFD_SP_clov[1][elem_ij0];
		lctmp0 += _OFFD_SP_clov[1][elem_ijb0];

		_SP_B[j0][i0] = lctmp0;
		_SP_B[i0][j0] = adj(lctmp0);
	      
		elem_ij0 = elem_ij0 + 1;
		elem_ijb0 = elem_ijb0 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;

  case 3:
    /*# gamma(  12)  -i  0  0  0            # ( 0011 )  --> 3 */
    /*#               0  i  0  0 */
    /*#               0  0 -i  0 */
    /*#               0  0  0  i */
    /*# From diagonal part */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp3;
	InnerReal lr_zero3;
	InnerReal lrtmp3;
	Complex sqm3;
	int i3;
	int j3;
	int elem_ij3;
	int elem_ijb3;
  
                          
	lr_zero3 = 0;
	FILL(sqm3, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i3 = 0; i3 < Nc; ++i3)
	    {
	      lrtmp3 = _DIAG_SP_clov[0][i3+Nc];
	      lrtmp3 -= _DIAG_SP_clov[0][i3];
	      lrtmp3 -= _DIAG_SP_clov[1][i3];
	      lrtmp3 += _DIAG_SP_clov[1][i3+Nc];
	      _SP_B[i3][i3] = cmplx(lr_zero3,lrtmp3);
	    }
	
	    /*# From lower triangular portion */
	    elem_ij3 = 0;
	    for(i3 = 1; i3 < Nc; ++i3)
	    {
	      elem_ijb3 = (i3+Nc)*(i3+Nc-1)/2 + Nc;
	      for(j3 = 0; j3 < i3; ++j3)
	      {
		lctmp3 = _OFFD_SP_clov[0][elem_ijb3];
		lctmp3 -= _OFFD_SP_clov[0][elem_ij3];
		lctmp3 -= _OFFD_SP_clov[1][elem_ij3];
		lctmp3 += _OFFD_SP_clov[1][elem_ijb3];

		_SP_B[j3][i3] = lctmp3 * sqm3;
		_SP_B[i3][j3] = adj(lctmp3) * sqm3;
	    
		elem_ij3 = elem_ij3 + 1;
		elem_ijb3 = elem_ijb3 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;

  case 5:
    /*# gamma(  13)   0 -1  0  0            # ( 0101 )  --> 5 */
    /*#               1  0  0  0 */
    /*#               0  0  0 -1 */
    /*#               0  0  1  0 */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp5;
	InnerReal lr_zero5;
	InnerReal lrtmp5;
	Complex sqm5;
	int i5;
	int j5;
	int elem_ij5;
	int elem_ji5;
  
                          
	lr_zero5 = 0;
	FILL(sqm5, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i5 = 0; i5 < Nc; ++i5)
	    {
	      elem_ij5 = (i5+Nc)*(i5+Nc-1)/2;
	      for(j5 = 0; j5 < Nc; ++j5)
	      {
		elem_ji5 = (j5+Nc)*(j5+Nc-1)/2 + i5;

		lctmp5 = adj(_OFFD_SP_clov[0][elem_ji5]);
		lctmp5 -= _OFFD_SP_clov[0][elem_ij5];
		lctmp5 += adj(_OFFD_SP_clov[1][elem_ji5]);
		lctmp5 -= _OFFD_SP_clov[1][elem_ij5];

		_SP_B[j5][i5] = lctmp5;

		elem_ij5 = elem_ij5 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;

  case 6:
    /*# gamma(  23)   0 -i  0  0            # ( 0110 )  --> 6 */
    /*#              -i  0  0  0 */
    /*#               0  0  0 -i */
    /*#               0  0 -i  0 */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp6;
	InnerReal lr_zero6;
	InnerReal lrtmp6;
	Complex sqm6;
	int i6;
	int j6;
	int elem_ij6;
	int elem_ji6;
  
                          
	lr_zero6 = 0;
	FILL(sqm6, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i6 = 0; i6 < Nc; ++i6)
	    {
	      elem_ij6 = (i6+Nc)*(i6+Nc-1)/2;
	      for(j6 = 0; j6 < Nc; ++j6)
	      {
		elem_ji6 = (j6+Nc)*(j6+Nc-1)/2 + i6;

		lctmp6 = adj(_OFFD_SP_clov[0][elem_ji6]);
		lctmp6 += _OFFD_SP_clov[0][elem_ij6];
		lctmp6 += adj(_OFFD_SP_clov[1][elem_ji6]);
		lctmp6 += _OFFD_SP_clov[1][elem_ij6];

		_SP_B[j6][i6] = -(lctmp6 * sqm6);

		elem_ij6 = elem_ij6 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;

  case 9:
    /*# gamma(  14)   0  i  0  0            # ( 1001 )  --> 9 */
    /*#               i  0  0  0 */
    /*#               0  0  0 -i */
    /*#               0  0 -i  0 */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp9;
	InnerReal lr_zero9;
	InnerReal lrtmp9;
	Complex sqm9;
	int i9;
	int j9;
	int elem_ij9;
	int elem_ji9;
  
                          
	lr_zero9 = 0;
	FILL(sqm9, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i9 = 0; i9 < Nc; ++i9)
	    {
	      elem_ij9 = (i9+Nc)*(i9+Nc-1)/2;
	      for(j9 = 0; j9 < Nc; ++j9)
	      {
		elem_ji9 = (j9+Nc)*(j9+Nc-1)/2 + i9;

		lctmp9 = adj(_OFFD_SP_clov[0][elem_ji9]);
		lctmp9 += _OFFD_SP_clov[0][elem_ij9];
		lctmp9 -= adj(_OFFD_SP_clov[1][elem_ji9]);
		lctmp9 -= _OFFD_SP_clov[1][elem_ij9];

		_SP_B[j9][i9] = lctmp9 * sqm9;

		elem_ij9 = elem_ij9 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;

  case 10:
    /*# gamma(  24)   0 -1  0  0            # ( 1010 )  --> 10 */
    /*#               1  0  0  0 */
    /*#               0  0  0  1 */
    /*#               0  0 -1  0 */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp10;
	InnerReal lr_zero10;
	InnerReal lrtmp10;
	Complex sqm10;
	int i10;
	int j10;
	int elem_ij10;
	int elem_ji10;
  
                          
	lr_zero10 = 0;
	FILL(sqm10, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i10 = 0; i10 < Nc; ++i10)
	    {
	      elem_ij10 = (i10+Nc)*(i10+Nc-1)/2;
	      for(j10 = 0; j10 < Nc; ++j10)
	      {
		elem_ji10 = (j10+Nc)*(j10+Nc-1)/2 + i10;

		lctmp10 = adj(_OFFD_SP_clov[0][elem_ji10]);
		lctmp10 -= _OFFD_SP_clov[0][elem_ij10];
		lctmp10 -= adj(_OFFD_SP_clov[1][elem_ji10]);
		lctmp10 += _OFFD_SP_clov[1][elem_ij10];

		_SP_B[j10][i10] = lctmp10;

		elem_ij10 = elem_ij10 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;
    
  case 12:
    /*# gamma(  34)   i  0  0  0            # ( 1100 )  --> 12 */
    /*#               0 -i  0  0 */
    /*#               0  0 -i  0 */
    /*#               0  0  0  i */
    /*# From diagonal part */
    SPECIAL(clov,B)
      {
	InnerComplex lctmp12;
	InnerReal lr_zero12;
	InnerReal lrtmp12;
	Complex sqm12;
	int i12;
	int j12;
	int elem_ij12;
	int elem_ijb12;
  
                          
	lr_zero12 = 0;
	FILL(sqm12, I);
  
	SPECIAL_LOOP_START()
	  {    
	    for(i12 = 0; i12 < Nc; ++i12)
	    {
	      lrtmp12 = _DIAG_SP_clov[0][i12];
	      lrtmp12 -= _DIAG_SP_clov[0][i12+Nc];
	      lrtmp12 -= _DIAG_SP_clov[1][i12];
	      lrtmp12 += _DIAG_SP_clov[1][i12+Nc];
	      _SP_B[i12][i12] = cmplx(lr_zero12,lrtmp12);
	    }
    
	    /*# From lower triangular portion */
	    elem_ij12 = 0;
	    for(i12 = 1; i12 < Nc; ++i12)
	    {
	      elem_ijb12 = (i12+Nc)*(i12+Nc-1)/2 + Nc;
	      for(j12 = 0; j12 < i12; ++j12)
	      {
		lctmp12 = _OFFD_SP_clov[0][elem_ij12];
		lctmp12 -= _OFFD_SP_clov[0][elem_ijb12];
		lctmp12 -= _OFFD_SP_clov[1][elem_ij12];
		lctmp12 += _OFFD_SP_clov[1][elem_ijb12];
	
		_SP_B[j12][i12] = lctmp12 * sqm12;
		_SP_B[i12][j12] = adj(lctmp12) * sqm12;
	
		elem_ij12 = elem_ij12 + 1;
		elem_ijb12 = elem_ijb12 + 1;
	      }
	    }
	  }
	SPECIAL_LOOP_END()

	  }
    SPECIAL_END(clov,B)
      break;
    
  default:
    B = 0;
  }
  

  END_CODE("subroutine");;
}
