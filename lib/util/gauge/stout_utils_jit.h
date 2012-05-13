

      void getFsAndBsJIT( GetFsAndBsArgs* arg )
      {
	const LatticeColorMatrix& Q = arg->Q;
	const LatticeColorMatrix& QQ = arg->QQ;
	multi1d<LatticeComplex>& f = arg->f;
	multi1d<LatticeComplex>& b1 = arg->b1;
	multi1d<LatticeComplex>& b2 = arg->b2;
	bool dobs=arg->dobs;

        while (1) {

	  static QDPJit<>::SharedLibEntry sharedLibEntry;
	  static MapVolumes*              mapVolumes;
	  static string                   strId;
	  static string                   prg;

	  QDPJitArgs cudaArgs;

	  string strREALT;
	  string codeQ,codeQQ;

          string codeF[3];
          string codeB1[3];
          string codeB2[3];


	  getTypeString( strREALT , REAL(0) );

	  if (!getCodeString( codeQ , Q , "site", cudaArgs )) { QDP_info("getFsAndBsJIT: could not cache Q"); break;  }      
	  if (!getCodeString( codeQQ , QQ , "site", cudaArgs )) { QDP_info("getFsAndBsJIT: could not cache QQ"); break;  }
	  if (!getCodeString( codeF[0] , f[0] , "site", cudaArgs )) { QDP_info("getFsAndBsJIT: could not cache f[0]"); break;}      
	  if (!getCodeString( codeF[1] , f[1] , "site", cudaArgs )) { QDP_info("getFsAndBsJIT: could not cache f[2]"); break;}
	  if (!getCodeString( codeF[2] , f[2] , "site", cudaArgs )) { QDP_info("getFsAndBsJIT: could not cache f[3]"); break;}
	  if (dobs) {
	    if (!getCodeString( codeB1[0] , b1[0] , "site", cudaArgs )) { QDP_info("getFsAndBsJIT: could not cache b1[0]"); break;} 
	    if (!getCodeString( codeB1[1] , b1[1] , "site", cudaArgs )) { QDP_info("getFsAndBsJIT: could not cache b1[2]"); break;}
	    if (!getCodeString( codeB1[2] , b1[2] , "site", cudaArgs )) { QDP_info("getFsAndBsJIT: could not cache b1[3]"); break;}
	    if (!getCodeString( codeB2[0] , b2[0] , "site", cudaArgs )) { QDP_info("getFsAndBsJIT: could not cache b2[0]"); break;} 
	    if (!getCodeString( codeB2[1] , b2[1] , "site", cudaArgs )) { QDP_info("getFsAndBsJIT: could not cache b2[2]"); break;}
	    if (!getCodeString( codeB2[2] , b2[2] , "site", cudaArgs )) { QDP_info("getFsAndBsJIT: could not cache b2[3]"); break;}
	  }

	  if (!mapVolumes) {
	    ostringstream osId;
	    osId << "getFsAndBsJIT1 "  << strREALT << " dobs=" << dobs;
	    strId = osId.str();
	    xmlready(strId);
#ifdef GPU_DEBUG_DEEP
	    cout << "strId = " << strId << endl;
#endif

	    std::ostringstream sprg;

	    sprg << " typedef " << strREALT << " REAL;" << endl;
	    sprg << " const int Nc = " << Nc << ";" << endl;

	    sprg << "int site = blockDim.x * blockIdx.x + blockDim.x * gridDim.x * blockIdx.y + threadIdx.x;" << endl;
	    sprg << "if (site >= args->numSiteTable) return;" << endl;

	    sprg << "   // Get the traces " << endl;
	    sprg << "   PColorMatrix<QDP::RComplex<REAL>, Nc>  Q_site = " << codeQ << ".elem(); " << endl;
	    sprg << "   PColorMatrix<QDP::RComplex<REAL>, Nc>  QQ_site = " << codeQQ << ".elem(); " << endl;
	    sprg << "   PColorMatrix<QDP::RComplex<REAL>, Nc>  QQQ = QQ_site*Q_site; " << endl;

	    sprg << " 	 " << endl;
	    sprg << "   PScalar<PScalar<RScalar<REAL> > > trQQQ;  " << endl;
	    sprg << "   trQQQ  = realTrace(QQQ); " << endl;
	    sprg << "   PScalar<PScalar<RScalar<REAL> > > trQQ; " << endl;
	    sprg << "   trQQ   = realTrace(QQ_site); " << endl;
	    sprg << " 	 " << endl;
	    sprg << "   REAL c0    = ((REAL)1/(REAL)3) * trQQQ.elem().elem().elem();  // eq 13 " << endl;
	    sprg << "   REAL c1    = ((REAL)1/(REAL)2) * trQQ.elem().elem().elem();	 // eq 15  " << endl;
	    sprg << " 	 " << endl;
	    sprg << " 	 " << endl;
	    sprg << "   if( c1 < 4.0e-3  )  " << endl;
	    sprg << "     { // RGE: set to 4.0e-3 (CM uses this value). I ran into nans with 1.0e-4 " << endl;
	    sprg << "       // ================================================================================ " << endl;
	    sprg << "       //  " << endl;
	    sprg << "       // Corner Case 1: if c1 < 1.0e-4 this implies c0max ~ 3x10^-7 " << endl;
	    sprg << "       //    and in this case the division c0/c0max in arccos c0/c0max can be undefined " << endl;
	    sprg << "       //    and produce NaN's " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       // In this case what we can do is get the f-s a different way. We go back to basics: " << endl;
	    sprg << "       // " << endl;
	    sprg << "       // We solve (using maple) the matrix equations using the eigenvalues  " << endl;
	    sprg << "       // " << endl;
	    sprg << "       //  [ 1, q_1, q_1^2 ] [ f_0 ]       [ exp( iq_1 ) ] " << endl;
	    sprg << "       //  [ 1, q_2, q_2^2 ] [ f_1 ]   =   [ exp( iq_2 ) ] " << endl;
	    sprg << "       //  [ 1, q_3, q_3^2 ] [ f_2 ]       [ exp( iq_3 ) ] " << endl;
	    sprg << "       // " << endl;
	    sprg << "       // with q_1 = 2 u w, q_2 = -u + w, q_3 = - u - w " << endl;
	    sprg << "       //  " << endl;
	    sprg << "       // with u and w defined as  u = sqrt( c_1/ 3 ) cos (theta/3) " << endl;
	    sprg << "       //                     and  w = sqrt( c_1 ) sin (theta/3) " << endl;
	    sprg << "       //                          theta = arccos ( c0 / c0max ) " << endl;
	    sprg << "       // leaving c0max as a symbol. " << endl;
	    sprg << "       // " << endl;
	    sprg << "       //  we then expand the resulting f_i as a series around c0 = 0 and c1 = 0 " << endl;
	    sprg << "       //  and then substitute in c0max = 2 ( c_1/ 3)^(3/2) " << endl;
	    sprg << "       //   " << endl;
	    sprg << "       //  we then convert the results to polynomials and take the real and imaginary parts: " << endl;
	    sprg << "       //  we get at the end of the day (to low order) " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       //                  1    2  " << endl;
	    sprg << "       //   f0[re] := 1 - --- c0  + h.o.t " << endl;
	    sprg << "       //                 720      " << endl;
	    sprg << "       // " << endl;
	    sprg << "       //	         1       1           1        2  " << endl;
	    sprg << "       //   f0[im] := - - c0 + --- c0 c1 - ---- c0 c1   + h.o.t " << endl;
	    sprg << "       //               6      120         5040         " << endl;
	    sprg << "       // " << endl;
	    sprg << "       // " << endl;
	    sprg << "       //             1        1            1        2  " << endl;
	    sprg << "       //   f1[re] := -- c0 - --- c0 c1 + ----- c0 c1  +  h.o.t " << endl;
	    sprg << "       //             24      360         13440        f " << endl;
	    sprg << "       // " << endl;
	    sprg << "       //                 1       1    2    1     3    1     2 " << endl;
	    sprg << "       //   f1[im] := 1 - - c1 + --- c1  - ---- c1  - ---- c0   + h.o.t " << endl;
	    sprg << "       //                 6      120       5040       5040 " << endl;
	    sprg << "       // " << endl;
	    sprg << "       //               1   1        1    2     1     3     1     2 " << endl;
	    sprg << "       //   f2[re] := - - + -- c1 - --- c1  + ----- c1  + ----- c0  + h.o.t " << endl;
	    sprg << "       //               2   24      720       40320       40320     " << endl;
	    sprg << "       // " << endl;
	    sprg << "       //              1        1              1        2 " << endl;
	    sprg << "       //   f2[im] := --- c0 - ---- c0 c1 + ------ c0 c1  + h.o.t " << endl;
	    sprg << "       //             120      2520         120960 " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       //  We then express these using Horner's rule for more stable evaluation. " << endl;
	    sprg << "       //  " << endl;
	    sprg << "       //  to get the b-s we use the fact that " << endl;
	    sprg << "       //                                      b2_i = d f_i / d c0 " << endl;
	    sprg << "       //                                 and  b1_i = d f_i / d c1 " << endl;
	    sprg << "       // " << endl;
	    sprg << "       //  where the derivatives are partial derivativs " << endl;
	    sprg << "       // " << endl;
	    sprg << "       //  And we just differentiate the polynomials above (keeping the same level " << endl;
	    sprg << "       //  of truncation) and reexpress that as Horner's rule " << endl;
	    sprg << "       //  " << endl;
	    sprg << "       //  This clearly also handles the case of a unit gauge as no c1, u etc appears in the  " << endl;
	    sprg << "       //  denominator and the arccos is never taken. In this case, we have the results in  " << endl;
	    sprg << "       //  the raw c0, c1 form and we don't need to flip signs and take complex conjugates. " << endl;
	    sprg << "       // " << endl;
	    sprg << "       //  I checked the expressions below by taking the difference between the Horner forms " << endl;
	    sprg << "       //  below from the expanded forms (and their derivatives) above and checking for the " << endl;
	    sprg << "       //  differences to be zero. At this point in time maple seems happy. " << endl;
	    sprg << "       //  ================================================================================== " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       " << codeF[0] << ".elem().elem().real() = 1.0-c0*c0/720.0; " << endl;
	    sprg << "       " << codeF[0] << ".elem().elem().imag() =  -(c0/6.0)*(1.0-(c1/20.0)*(1.0-(c1/42.0))) ; " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       " << codeF[1] << ".elem().elem().real() =  c0/24.0*(1.0-c1/15.0*(1.0-3.0*c1/112.0)) ; " << endl;
	    sprg << "       " << codeF[1] << ".elem().elem().imag() =  1.0-c1/6.0*(1.0-c1/20.0*(1.0-c1/42.0))-c0*c0/5040.0 ; " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       " << codeF[2] << ".elem().elem().real() = 0.5*(-1.0+c1/12.0*(1.0-c1/30.0*(1.0-c1/56.0))+c0*c0/20160.0); " << endl;
	    sprg << "       " << codeF[2] << ".elem().elem().imag() = 0.5*(c0/60.0*(1.0-c1/21.0*(1.0-c1/48.0))); " << endl;
	    sprg << " 	   " << endl;
	    if (dobs) {
	      sprg << " 	//  partial f0/ partial c0 " << endl;
	      sprg << " 	" << codeB2[0] << ".elem().elem().real() = -c0/360.0; " << endl;
	      sprg << " 	" << codeB2[0] << ".elem().elem().imag() =  -(1.0/6.0)*(1.0-(c1/20.0)*(1.0-c1/42.0)); " << endl;
	      sprg << " 	     " << endl;
	      sprg << " 	// partial f0 / partial c1 " << endl;
	      sprg << " 	// " << endl;
	      sprg << " 	" << codeB1[0] << ".elem().elem().real() = 0; " << endl;
	      sprg << " 	" << codeB1[0] << ".elem().elem().imag() = (c0/120.0)*(1.0-c1/21.0); " << endl;
	      sprg << " 	     " << endl;
	      sprg << " 	// partial f1 / partial c0 " << endl;
	      sprg << " 	// " << endl;
	      sprg << " 	" << codeB2[1] << ".elem().elem().real() = (1.0/24.0)*(1.0-c1/15.0*(1.0-3.0*c1/112.0)); " << endl;
	      sprg << " 	" << codeB2[1] << ".elem().elem().imag() = -c0/2520.0; " << endl;
	      sprg << " 	     " << endl;
	      sprg << " 	     " << endl;
	      sprg << " 	// partial f1 / partial c1 " << endl;
	      sprg << " 	" << codeB1[1] << ".elem().elem().real() = -c0/360.0*(1.0 - 3.0*c1/56.0 ); " << endl;
	      sprg << " 	" << codeB1[1] << ".elem().elem().imag() = -1.0/6.0*(1.0-c1/10.0*(1.0-c1/28.0)); " << endl;
	      sprg << " 	     " << endl;
	      sprg << " 	// partial f2/ partial c0 " << endl;
	      sprg << " 	" << codeB2[2] << ".elem().elem().real() = 0.5*c0/10080.0; " << endl;
	      sprg << " 	" << codeB2[2] << ".elem().elem().imag() = 0.5*(  1.0/60.0*(1.0-c1/21.0*(1.0-c1/48.0)) ); " << endl;
	      sprg << " 	     " << endl;
	      sprg << " 	// partial f2/ partial c1 " << endl;
	      sprg << " 	" << codeB1[2] << ".elem().elem().real() = 0.5*(  1.0/12.0*(1.0-(2.0*c1/30.0)*(1.0-3.0*c1/112.0)) );  " << endl;
	      sprg << " 	" << codeB1[2] << ".elem().elem().imag() = 0.5*( -c0/1260.0*(1.0-c1/24.0) ); " << endl;
	      sprg << " 	     " << endl;
	    }
	    sprg << "     } " << endl;
	    sprg << "   else  " << endl;
	    sprg << "     {  " << endl;
	    sprg << "       // =================================================================================== " << endl;
	    sprg << "       // Normal case: Do as per paper " << endl;
	    sprg << "       // =================================================================================== " << endl;
	    sprg << "       bool c0_negativeP = c0 < 0; " << endl;
	    sprg << "       REAL c0abs = fabs((double)c0); " << endl;
	    sprg << "       REAL c0max = 2*pow( (double)(c1/(double)3), (double)1.5); " << endl;
	    sprg << "       REAL theta; " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       // ====================================================================================== " << endl;
	    sprg << "       // Now work out theta. In the paper the case where c0 -> c0max even when c1 is reasonable  " << endl;
	    sprg << "       // Has never been considered, even though it can arise and can cause the arccos function " << endl;
	    sprg << "       // to fail " << endl;
	    sprg << "       // Here we handle it with series expansion " << endl;
	    sprg << "       // ===================================================================================== " << endl;
	    sprg << "       REAL eps = (c0max - c0abs)/c0max; " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       if( eps < 0 ) { " << endl;
	    sprg << " 	// =============================================================================== " << endl;
	    sprg << " 	// Corner Case 2: Handle case when c0abs is bigger than c0max.  " << endl;
	    sprg << " 	// This can happen only when there is a rounding error in the ratio, and that the  " << endl;
	    sprg << " 	// ratio is really 1. This implies theta = 0 which we'll just set. " << endl;
	    sprg << " 	// =============================================================================== " << endl;
	    sprg << " 	theta = 0; " << endl;
	    sprg << "       } " << endl;
	    sprg << "       else if ( eps < 1.0e-3 ) { " << endl;
	    sprg << " 	// =============================================================================== " << endl;
	    sprg << " 	// Corner Case 3: c0->c0max even though c1 may be actually quite reasonable. " << endl;
	    sprg << " 	// The ratio |c0|/c0max -> 1 but is still less than one, so that a  " << endl;
	    sprg << " 	// series expansion is possible. " << endl;
	    sprg << " 	// SERIES of acos(1-epsilon): Good to O(eps^6) or with this cutoff to O(10^{-18}) Computed with Maple. " << endl;
	    sprg << " 	//  BTW: 1-epsilon = 1 - (c0max-c0abs)/c0max = 1-(1 - c0abs/c0max) = +c0abs/c0max " << endl;
	    sprg << " 	// " << endl;
	    sprg << " 	// =============================================================================== " << endl;
	    sprg << " 	REAL sqtwo = sqrt((REAL)2); " << endl;
	    sprg << " 	     " << endl;
	    sprg << " 	theta = sqtwo*sqrt(eps)*( 1.0 + ( (1/(REAL)12) + ( (3/(REAL)160) + ( (5/(REAL)896) + ( (35/(REAL)18432) + (63/(REAL)90112)*eps ) *eps) *eps) *eps) *eps); " << endl;
	    sprg << " 	     " << endl;
	    sprg << "       }  " << endl;
	    sprg << "       else {   " << endl;
	    sprg << " 	//  " << endl;
	    sprg << " 	theta = acos( c0abs/c0max ); " << endl;
	    sprg << "       } " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       REAL f_site_re[3]; " << endl;
	    sprg << "       REAL f_site_im[3]; " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       REAL b1_site_re[3]; " << endl;
	    sprg << "       REAL b1_site_im[3]; " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       REAL b2_site_re[3]; " << endl;
	    sprg << "       REAL b2_site_im[3]; " << endl;
	    sprg << " 	   " << endl;
	    sprg << " 	   " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       REAL u = sqrt(c1/3)*cos(theta/3); " << endl;
	    sprg << "       REAL w = sqrt(c1)*sin(theta/3); " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       REAL u_sq = u*u; " << endl;
	    sprg << "       REAL w_sq = w*w; " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       REAL xi0,xi1; " << endl;
	    sprg << "       { " << endl;
	    sprg << " 	bool w_smallP  = fabs(w) < 0.05; " << endl;
	    sprg << " 	if( w_smallP ) {  " << endl;
	    sprg << " 	  xi0 = (REAL)1 - ((REAL)1/(REAL)6)*w_sq*( 1 - ((REAL)1/(REAL)20)*w_sq*( (REAL)1 - ((REAL)1/(REAL)42)*w_sq ) ); " << endl;
	    sprg << " 	} " << endl;
	    sprg << " 	else { " << endl;
	    sprg << " 	  xi0 = sin(w)/w; " << endl;
	    sprg << " 	} " << endl;
	    sprg << " 	     " << endl;
	    if (dobs) {
	      sprg << " 	  if( w_smallP  ) {  " << endl;
	      sprg << " 	    xi1 = -1*( ((REAL)1/(REAL)3) - ((REAL)1/(REAL)30)*w_sq*( (REAL)1 - ((REAL)1/(REAL)28)*w_sq*( (REAL)1 - ((REAL)1/(REAL)54)*w_sq ) ) ); " << endl;
	      sprg << " 	  } " << endl;
	      sprg << " 	  else {  " << endl;
	      sprg << " 	    xi1 = cos(w)/w_sq - sin(w)/(w_sq*w); " << endl;
	      sprg << " 	  } " << endl;
	    }
	    sprg << "       } " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       REAL cosu = cos(u); " << endl;
	    sprg << "       REAL sinu = sin(u); " << endl;
	    sprg << "       REAL cosw = cos(w); " << endl;
	    sprg << "       REAL sinw = sin(w); " << endl;
	    sprg << "       REAL sin2u = sin(2*u); " << endl;
	    sprg << "       REAL cos2u = cos(2*u); " << endl;
	    sprg << "       REAL ucosu = u*cosu; " << endl;
	    sprg << "       REAL usinu = u*sinu; " << endl;
	    sprg << "       REAL ucos2u = u*cos2u; " << endl;
	    sprg << "       REAL usin2u = u*sin2u; " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       REAL denum = (REAL)9*u_sq - w_sq; " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       { " << endl;
	    sprg << " 	REAL subexp1 = u_sq - w_sq; " << endl;
	    sprg << " 	REAL subexp2 = 8*u_sq*cosw; " << endl;
	    sprg << " 	REAL subexp3 = (3*u_sq + w_sq)*xi0; " << endl;
	    sprg << " 	     " << endl;
	    sprg << " 	f_site_re[0] = ( (subexp1)*cos2u + cosu*subexp2 + 2*usinu*subexp3 ) / denum ; " << endl;
	    sprg << " 	f_site_im[0] = ( (subexp1)*sin2u - sinu*subexp2 + 2*ucosu*subexp3 ) / denum ; " << endl;
	    sprg << "       } " << endl;
	    sprg << "       { " << endl;
	    sprg << " 	REAL subexp = (3*u_sq -w_sq)*xi0; " << endl;
	    sprg << " 	     " << endl;
	    sprg << " 	f_site_re[1] = (2*(ucos2u - ucosu*cosw)+subexp*sinu)/denum; " << endl;
	    sprg << " 	f_site_im[1] = (2*(usin2u + usinu*cosw)+subexp*cosu)/denum; " << endl;
	    sprg << "       } " << endl;
	    sprg << " 	   " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       { " << endl;
	    sprg << " 	REAL subexp=3*xi0; " << endl;
	    sprg << " 	     " << endl;
	    sprg << " 	f_site_re[2] = (cos2u - cosu*cosw -usinu*subexp) /denum ; " << endl;
	    sprg << " 	f_site_im[2] = (sin2u + sinu*cosw -ucosu*subexp) /denum ; " << endl;
	    sprg << "       } " << endl;
	    sprg << " 	   " << endl;
	    if (dobs) {
	      sprg << " 	{ " << endl;
	      sprg << " 	  REAL r_1_re[3]; " << endl;
	      sprg << " 	  REAL r_1_im[3]; " << endl;
	      sprg << " 	  REAL r_2_re[3]; " << endl;
	      sprg << " 	  REAL r_2_im[3]; " << endl;
	      sprg << " 	       " << endl;
	      sprg << " 	  //	  r_1[0]=Double(2)*cmplx(u, u_sq-w_sq)*exp2iu " << endl;
	      sprg << " 	  //          + 2.0*expmiu*( cmplx(8.0*u*cosw, -4.0*u_sq*cosw) " << endl;
	      sprg << " 	  //	      + cmplx(u*(3.0*u_sq+w_sq),9.0*u_sq+w_sq)*xi0 ); " << endl;
	      sprg << " 	  { " << endl;
	      sprg << " 	    REAL subexp1 = u_sq - w_sq; " << endl;
	      sprg << " 	    REAL subexp2 =  8*cosw + (3*u_sq + w_sq)*xi0 ; " << endl;
	      sprg << " 	    REAL subexp3 =  4*u_sq*cosw - (9*u_sq + w_sq)*xi0 ; " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	    r_1_re[0] = 2*(ucos2u - sin2u *(subexp1)+ucosu*( subexp2 )- sinu*( subexp3 ) ); " << endl;
	      sprg << " 	    r_1_im[0] = 2*(usin2u + cos2u *(subexp1)-usinu*( subexp2 )- cosu*( subexp3 ) ); " << endl;
	      sprg << " 	  } " << endl;
	      sprg << " 	       " << endl;
	      sprg << " 	  // r_1[1]=cmplx(2.0, 4.0*u)*exp2iu + expmiu*cmplx(-2.0*cosw-(w_sq-3.0*u_sq)*xi0,2.0*u*cosw+6.0*u*xi0); " << endl;
	      sprg << " 	  { " << endl;
	      sprg << " 	    REAL subexp1 = cosw+3*xi0; " << endl;
	      sprg << " 	    REAL subexp2 = 2*cosw + xi0*(w_sq - 3*u_sq); " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	    r_1_re[1] = 2*((cos2u - 2*usin2u) + usinu*( subexp1 )) - cosu*( subexp2 ); " << endl;
	      sprg << " 	    r_1_im[1] = 2*((sin2u + 2*ucos2u) + ucosu*( subexp1 )) + sinu*( subexp2 ); " << endl;
	      sprg << " 	  } " << endl;
	      sprg << " 	       " << endl;
	      sprg << " 	       " << endl;
	      sprg << " 	  // r_1[2]=2.0*timesI(exp2iu)  +expmiu*cmplx(-3.0*u*xi0, cosw-3*xi0); " << endl;
	      sprg << " 	  { " << endl;
	      sprg << " 	    REAL subexp = cosw - 3*xi0; " << endl;
	      sprg << " 	    r_1_re[2] = -2*sin2u -3*ucosu*xi0 + sinu*( subexp ); " << endl;
	      sprg << " 	    r_1_im[2] = 2*cos2u  +3*usinu*xi0 + cosu*( subexp ); " << endl;
	      sprg << " 	  } " << endl;
	      sprg << " 	       " << endl;
	      sprg << " 	       " << endl;
	      sprg << " 	  //r_2[0]=-2.0*exp2iu + 2*cmplx(0,u)*expmiu*cmplx(cosw+xi0+3*u_sq*xi1, " << endl;
	      sprg << " 	  //						 4*u*xi0); " << endl;
	      sprg << " 	  { " << endl;
	      sprg << " 	    REAL subexp = cosw + xi0 + 3*u_sq*xi1; " << endl;
	      sprg << " 	    r_2_re[0] = -2*(cos2u + u*( 4*ucosu*xi0 - sinu*(subexp )) ); " << endl;
	      sprg << " 	    r_2_im[0] = -2*(sin2u - u*( 4*usinu*xi0 + cosu*(subexp )) ); " << endl;
	      sprg << " 	  } " << endl;
	      sprg << " 	       " << endl;
	      sprg << " 	       " << endl;
	      sprg << " 	  // r_2[1]= expmiu*cmplx(cosw+xi0-3.0*u_sq*xi1, 2.0*u*xi0); " << endl;
	      sprg << " 	  // r_2[1] = timesMinusI(r_2[1]); " << endl;
	      sprg << " 	  { " << endl;
	      sprg << " 	    REAL subexp =  cosw + xi0 - 3*u_sq*xi1; " << endl;
	      sprg << " 	    r_2_re[1] =  2*ucosu*xi0 - sinu*( subexp ) ; " << endl;
	      sprg << " 	    r_2_im[1] = -2*usinu*xi0 - cosu*( subexp ) ; " << endl;
	      sprg << " 	  } " << endl;
	      sprg << " 	       " << endl;
	      sprg << " 	  //r_2[2]=expmiu*cmplx(xi0, -3.0*u*xi1); " << endl;
	      sprg << " 	  { " << endl;
	      sprg << " 	    REAL subexp = 3*xi1; " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	    r_2_re[2] =    cosu*xi0 - usinu*subexp ; " << endl;
	      sprg << " 	    r_2_im[2] = -( sinu*xi0 + ucosu*subexp ) ; " << endl;
	      sprg << " 	  }       " << endl;
	      sprg << " 	       " << endl;
	      sprg << " 	  REAL b_denum=2*denum*denum; " << endl;
	      sprg << " 	       " << endl;
	      sprg << " 	       " << endl;
	      sprg << " 	  for(int j=0; j < 3; j++) {  " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	    { " << endl;
	      sprg << " 	      REAL subexp1 = 2*u; " << endl;
	      sprg << " 	      REAL subexp2 = 3*u_sq - w_sq; " << endl;
	      sprg << " 	      REAL subexp3 = 2*(15*u_sq + w_sq); " << endl;
	      sprg << " 		   " << endl;
	      sprg << " 	      b1_site_re[j]=( subexp1*r_1_re[j] + subexp2*r_2_re[j] - subexp3*f_site_re[j] )/b_denum; " << endl;
	      sprg << " 	      b1_site_im[j]=( subexp1*r_1_im[j] + subexp2*r_2_im[j] - subexp3*f_site_im[j] )/b_denum; " << endl;
	      sprg << " 	    } " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	    {  " << endl;
	      sprg << " 	      REAL subexp1 = 3*u; " << endl;
	      sprg << " 	      REAL subexp2 = 24*u; " << endl;
	      sprg << " 		   " << endl;
	      sprg << " 	      b2_site_re[j]=( r_1_re[j]- subexp1*r_2_re[j] - subexp2 * f_site_re[j] )/b_denum; " << endl;
	      sprg << " 	      b2_site_im[j]=( r_1_im[j] -subexp1*r_2_im[j] - subexp2 * f_site_im[j] )/b_denum; " << endl;
	      sprg << " 	    } " << endl;
	      sprg << " 	  } " << endl;
	      sprg << "  " << endl;
	      sprg << " 	  // Now flip the coefficients of the b-s " << endl;
	      sprg << " 	  if( c0_negativeP )  " << endl;
	      sprg << " 	    { " << endl;
	      sprg << " 	      //b1_site[0] = conj(b1_site[0]); " << endl;
	      sprg << " 	      b1_site_im[0] *= -1; " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	      //b1_site[1] = -conj(b1_site[1]); " << endl;
	      sprg << " 	      b1_site_re[1] *= -1; " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	      //b1_site[2] = conj(b1_site[2]); " << endl;
	      sprg << " 	      b1_site_im[2] *= -1; " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	      //b2_site[0] = -conj(b2_site[0]); " << endl;
	      sprg << " 	      b2_site_re[0] *= -1; " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	      //b2_site[1] = conj(b2_site[1]); " << endl;
	      sprg << " 	      b2_site_im[1] *= -1; " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	      //b2_site[2] = -conj(b2_site[2]); " << endl;
	      sprg << " 	      b2_site_re[2] *= -1; " << endl;
	      sprg << " 	    } " << endl;
	      sprg << " 	       " << endl;
	      sprg << " 	  // Load back into the lattice sized object " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	    " << codeB1[0] << ".elem().elem().real() = b1_site_re[0]; " << endl;
	      sprg << " 	    " << codeB1[0] << ".elem().elem().imag() = b1_site_im[0]; " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	    " << codeB2[0] << ".elem().elem().real() = b2_site_re[0]; " << endl;
	      sprg << " 	    " << codeB2[0] << ".elem().elem().imag() = b2_site_im[0]; " << endl;

	      sprg << " 		 " << endl;
	      sprg << " 	    " << codeB1[1] << ".elem().elem().real() = b1_site_re[1]; " << endl;
	      sprg << " 	    " << codeB1[1] << ".elem().elem().imag() = b1_site_im[1]; " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	    " << codeB2[1] << ".elem().elem().real() = b2_site_re[1]; " << endl;
	      sprg << " 	    " << codeB2[1] << ".elem().elem().imag() = b2_site_im[1]; " << endl;

	      sprg << " 		 " << endl;
	      sprg << " 	    " << codeB1[2] << ".elem().elem().real() = b1_site_re[2]; " << endl;
	      sprg << " 	    " << codeB1[2] << ".elem().elem().imag() = b1_site_im[2]; " << endl;
	      sprg << " 		 " << endl;
	      sprg << " 	    " << codeB2[2] << ".elem().elem().real() = b2_site_re[2]; " << endl;
	      sprg << " 	    " << codeB2[2] << ".elem().elem().imag() = b2_site_im[2]; " << endl;

	      sprg << " 	   }" << endl;
	    }
	    sprg << " 	   " << endl;
	    sprg << "       // Now when everything is done flip signs of the b-s (can't do this before " << endl;
	    sprg << "       // as the unflipped f-s are needed to find the b-s " << endl;
	    sprg << " 	   " << endl;
	    sprg << "       if( c0_negativeP ) { " << endl;
	    sprg << " 	     " << endl;
	    sprg << " 	// f_site[0] = conj(f_site[0]); " << endl;
	    sprg << " 	f_site_im[0] *= -1; " << endl;
	    sprg << " 	     " << endl;
	    sprg << " 	//f_site[1] = -conj(f_site[1]); " << endl;
	    sprg << " 	f_site_re[1] *= -1; " << endl;
	    sprg << " 	     " << endl;
	    sprg << " 	//f_site[2] = conj(f_site[2]); " << endl;
	    sprg << " 	f_site_im[2] *= -1; " << endl;
	    sprg << " 	     " << endl;
	    sprg << "       } " << endl;
	    sprg << " 	 " << endl;
	    sprg << "       // Load back into the lattice sized object " << endl;
	    sprg << "      " << codeF[0] << ".elem().elem().real() = f_site_re[0];" << endl;
	    sprg << "      " << codeF[0] << ".elem().elem().imag() = f_site_im[0];" << endl;
	    sprg << "      " << codeF[1] << ".elem().elem().real() = f_site_re[1];" << endl;
	    sprg << "      " << codeF[1] << ".elem().elem().imag() = f_site_im[1];" << endl;
	    sprg << "      " << codeF[2] << ".elem().elem().real() = f_site_re[2];" << endl;
	    sprg << "      " << codeF[2] << ".elem().elem().imag() = f_site_im[2];" << endl;
	    sprg << " 	   " << endl;
	    sprg << "     } // End of if( corner_caseP ) else {} " << endl;




	    prg = sprg.str();

#ifdef GPU_DEBUG_DEEP
	    cout << "Cuda kernel code = " << endl << prg << endl << endl;
#endif
	  }

	  const int nodeSites = QDP::Layout::sitesOnNode();

	  cudaArgs.getArgs().numSiteTable = nodeSites;
	  cudaArgs.getArgs().indirection  = arg->dobs;
   
	  StopWatch watch0;
	  watch0.start();

	  if (!QDPJit<>::Instance()( strId , prg , cudaArgs , nodeSites , sharedLibEntry  , mapVolumes )) {
	    QDP_info("getFsAndBsJIT: call to cuda jitter failed");
	    break;
	  }

	  for ( int i = 0 ; i < 3 ; i++) {
	    f[i].setModified(true);
	    if (dobs) {
	      b1[i].setModified(true);
	      b2[i].setModified(true);
	    }
	  }
	      

	  watch0.stop();
	  DeviceStats::Instance().incMicroSecondsKernelExec( EVAL_LAT_LAT,watch0.getTimeInMicroseconds() );
	  DeviceStats::Instance().incEvalDev(EVAL_LAT_LAT);
	
	  QDP::theCacheDispatch::Instance().unLock();

	  return;
	}

	QDP::theCacheDispatch::Instance().unLock();
  
	QDP_error("getFsAndBsJIT: host not implemented");
      }



