#include "qdp.h"

#ifdef QDP_IS_QDPJIT

using namespace QDP;


CUfunction function_get_fs_bs_exec(CUfunction function, 
				   const LatticeColorMatrix& Q,
				   const LatticeColorMatrix& QQ,
				   multi1d<LatticeComplex>& f,
				   multi1d<LatticeComplex>& b1,
				   multi1d<LatticeComplex>& b2,
				   bool dobs)
{
  AddressLeaf addr_leaf;

  int junk_0 = forEach(Q, addr_leaf, NullCombine());
  int junk_1 = forEach(QQ, addr_leaf, NullCombine());
  int junk_2 = forEach(f[0], addr_leaf, NullCombine());
  int junk_3 = forEach(f[1], addr_leaf, NullCombine());
  int junk_4 = forEach(f[2], addr_leaf, NullCombine());
  int junk_5 = forEach(b1[0], addr_leaf, NullCombine());
  int junk_6 = forEach(b1[1], addr_leaf, NullCombine());
  int junk_7 = forEach(b1[2], addr_leaf, NullCombine());
  int junk_8 = forEach(b2[0], addr_leaf, NullCombine());
  int junk_9 = forEach(b2[1], addr_leaf, NullCombine());
  int junk_10= forEach(b2[2], addr_leaf, NullCombine());


  // lo <= idx < hi
  int lo = 0;
  int hi = Layout::sitesOnNode();
  unsigned short dobs_u8 = dobs ? 1 : 0;

  std::vector<void*> addr;

  addr.push_back( &lo );
  //std::cout << "addr lo = " << addr[0] << " lo=" << lo << "\n";

  addr.push_back( &hi );
  //std::cout << "addr hi = " << addr[1] << " hi=" << hi << "\n";

  addr.push_back( &dobs_u8 );

  int addr_dest=addr.size();
  for(int i=0; i < addr_leaf.addr.size(); ++i) {
    addr.push_back( &addr_leaf.addr[i] );
    //std::cout << "addr = " << addr_leaf.addr[i] << "\n";
  }

  jit_launch(function,hi-lo,addr);
}



WordREG<REAL> jit_constant( double f )
{
  return WordREG<REAL>(f);
}


CUfunction function_get_fs_bs_build(const LatticeColorMatrix& Q,
				    const LatticeColorMatrix& QQ,
				    multi1d<LatticeComplex>& f,
				    multi1d<LatticeComplex>& b1,
				    multi1d<LatticeComplex>& b2,
				    bool dobs)
{
  //std::cout << __PRETTY_FUNCTION__ << ": entering\n";

  CUfunction func;

  jit_start_new_function();

  jit_value r_lo     = jit_add_param( jit_ptx_type::s32 );
  jit_value r_hi     = jit_add_param( jit_ptx_type::s32 );
  jit_value r_dobs   = jit_add_param( jit_ptx_type::pred );
  jit_value r_nobs   = jit_ins_not( r_dobs );
      
  jit_value r_idx = jit_geom_get_linear_th_idx();  
      
  jit_value r_out_of_range       = jit_ins_ge( r_idx , r_hi );
  jit_ins_exit( r_out_of_range );
      
  ParamLeaf param_leaf( r_idx );
  
  typedef typename LeafFunctor<LatticeColorMatrix, ParamLeaf>::Type_t  LCMJIT;
  typedef typename LeafFunctor<LatticeComplex    , ParamLeaf>::Type_t  LCJIT;

  LCMJIT Q_jit(forEach(Q, param_leaf, TreeCombine()));
  LCMJIT QQ_jit(forEach(QQ, param_leaf, TreeCombine()));
  LCJIT  f0_jit(forEach(f[0], param_leaf, TreeCombine()));
  LCJIT  f1_jit(forEach(f[1], param_leaf, TreeCombine()));
  LCJIT  f2_jit(forEach(f[2], param_leaf, TreeCombine()));
  LCJIT  b10_jit(forEach(b1[0], param_leaf, TreeCombine()));
  LCJIT  b11_jit(forEach(b1[1], param_leaf, TreeCombine()));
  LCJIT  b12_jit(forEach(b1[2], param_leaf, TreeCombine()));
  LCJIT  b20_jit(forEach(b2[0], param_leaf, TreeCombine()));
  LCJIT  b21_jit(forEach(b2[1], param_leaf, TreeCombine()));
  LCJIT  b22_jit(forEach(b2[2], param_leaf, TreeCombine()));

  auto& Q_j  = Q_jit.elem(JitDeviceLayout::Coalesced);
  auto& QQ_j = QQ_jit.elem(JitDeviceLayout::Coalesced);

  auto& f0_j = f0_jit.elem(JitDeviceLayout::Coalesced);
  auto& f1_j = f1_jit.elem(JitDeviceLayout::Coalesced);
  auto& f2_j = f2_jit.elem(JitDeviceLayout::Coalesced);
  
  auto& b10_j = b10_jit.elem(JitDeviceLayout::Coalesced);
  auto& b11_j = b11_jit.elem(JitDeviceLayout::Coalesced);
  auto& b12_j = b12_jit.elem(JitDeviceLayout::Coalesced);
  auto& b20_j = b20_jit.elem(JitDeviceLayout::Coalesced);
  auto& b21_j = b21_jit.elem(JitDeviceLayout::Coalesced);
  auto& b22_j = b22_jit.elem(JitDeviceLayout::Coalesced);


  { 
    // Get the traces
    PColorMatrixREG< RComplexREG< WordREG<REAL> >, Nc>  Q_site = Q_j.elem();
    PColorMatrixREG< RComplexREG< WordREG<REAL> >, Nc>  QQ_site = QQ_j.elem();
    PColorMatrixREG< RComplexREG< WordREG<REAL> >, Nc>  QQQ = QQ_site*Q_site;
	
    // Real trQQQ; 
    // trQQQ.elem()  = realTrace(QQQ);
    // Real trQQ;
    // trQQ.elem()   = realTrace(QQ_site);

    // auto trQQQ = realTrace(QQQ);
    // auto trQQ  = realTrace(QQ_site);

    PScalarREG< RScalarREG< WordREG<REAL> > > trQQQ = realTrace(QQQ);
    PScalarREG< RScalarREG< WordREG<REAL> > > trQQ  = realTrace(QQ_site);

    WordREG<REAL> c0    = jit_constant((REAL)1/(REAL)3) * trQQQ.elem().elem();  // eq 13
    WordREG<REAL> c1    = jit_constant((REAL)1/(REAL)2) * trQQ.elem().elem();	 // eq 15 

    jit_label_t not_c1_lt;
    jit_label_t label_exit;
    jit_ins_branch( not_c1_lt , jit_ins_ge( c1.get_val() , jit_value( 4.0e-3 ) ) );

    { //    if( c1 < 4.0e-3  ) 
      f0_j.elem().elem().real() = jit_constant(1.0) - c0 * c0 / jit_constant(720.0);
      f0_j.elem().elem().imag() =  -( c0 / jit_constant(6.0) )*( jit_constant(1.0) -(c1/jit_constant(20.0))*(jit_constant(1.0)-(c1/jit_constant(42.0)))) ;

      f1_j.elem().elem().real() =  c0/jit_constant(24.0)*(jit_constant(1.0)-c1/jit_constant(15.0)*(jit_constant(1.0)-jit_constant(3.0)*c1/jit_constant(112.0))) ;
      f1_j.elem().elem().imag() =  jit_constant(1.0)-c1/jit_constant(6.0)*(jit_constant(1.0)-c1/jit_constant(20.0)*(jit_constant(1.0)-c1/jit_constant(42.0)))-c0*c0/jit_constant(5040.0);
	  
      f2_j.elem().elem().real() = jit_constant(0.5)*(jit_constant(-1.0)+c1/jit_constant(12.0)*(jit_constant(1.0)-c1/jit_constant(30.0)*(jit_constant(1.0)-c1/jit_constant(56.0)))+c0*c0/jit_constant(20160.0));
      f2_j.elem().elem().imag() = jit_constant(0.5)*(c0/jit_constant(60.0)*(jit_constant(1.0)-c1/jit_constant(21.0)*(jit_constant(1.0)-c1/jit_constant(48.0))));


      jit_label_t cont_0;
      jit_ins_branch( cont_0 , r_nobs );

      { // dobs

	b20_j.elem().elem().real() = -c0/jit_constant(360.0);
	b20_j.elem().elem().imag() =  -jit_constant(1.0/6.0)*(jit_constant(1.0)-(c1/jit_constant(20.0))*(jit_constant(1.0)-c1/jit_constant(42.0)));
	    
	// partial f0 / partial c1
	//
	b10_j.elem().elem().real() = jit_constant(0);
	b10_j.elem().elem().imag() = (c0/jit_constant(120.0))*(jit_constant(1.0)-c1/jit_constant(21.0));
	    
	// partial f1 / partial c0
	//
	b21_j.elem().elem().real() = jit_constant(1.0/24.0)*(jit_constant(1.0)-c1/jit_constant(15.0)*(jit_constant(1.0)-jit_constant(3.0)*c1/jit_constant(112.0)));
	b21_j.elem().elem().imag() = -c0/jit_constant(2520.0);
	    
	    
	// partial f1 / partial c1
	b11_j.elem().elem().real() = -c0/jit_constant(360.0)*(jit_constant(1.0) - jit_constant(3.0)*c1/jit_constant(56.0) );
	b11_j.elem().elem().imag() = -jit_constant(1.0/6.0)*(jit_constant(1.0)-c1/jit_constant(10.0)*(jit_constant(1.0)-c1/jit_constant(28.0)));
	    
	// partial f2/ partial c0
	b22_j.elem().elem().real() = jit_constant(0.5)*c0/jit_constant(10080.0);
	b22_j.elem().elem().imag() = jit_constant(0.5)*(  jit_constant(1.0/60.0)*(jit_constant(1.0)-c1/jit_constant(21.0)*(jit_constant(1.0)-c1/jit_constant(48.0))) );
	    
	// partial f2/ partial c1
	b12_j.elem().elem().real() = jit_constant(0.5)*(  jit_constant(1.0/12.0)*(jit_constant(1.0)-(jit_constant(2.0)*c1/jit_constant(30.0))*(jit_constant(1.0)-jit_constant(3.0)*c1/jit_constant(112.0))) ); 
	b12_j.elem().elem().imag() = jit_constant(0.5)*( -c0/jit_constant(1260.0)*(jit_constant(1.0)-c1/jit_constant(24.0)) );
      } // Dobs==true
      jit_ins_label( cont_0 );
      jit_ins_branch( label_exit );
    } // if (c1 < 4.0e-3 )
    jit_ins_label( not_c1_lt );


    jit_value c0_negativeP = jit_ins_lt( c0.get_val() , jit_value(0.0) );
    WordREG<REAL> c0abs = fabs(c0);
    WordREG<REAL> c0max = jit_constant(2.0) * pow( c1 / jit_constant(3.0) , jit_constant(1.5) );
    WordREG<REAL> theta;
    WordREG<REAL> eps = (c0max - c0abs)/c0max;


    jit_label_t cont_1;
    jit_label_t label_theta_exit;
    jit_ins_branch( cont_1 , jit_ins_ge( eps.get_val() , jit_value( 0.0 ) ) );
    // if( eps < 0 ) {
    // ===============================================================================
    // Corner Case 2: Handle case when c0abs is bigger than c0max. 
    // This can happen only when there is a rounding error in the ratio, and that the 
    // ratio is really 1. This implies theta = 0 which we'll just set.
    // ===============================================================================
    theta = jit_constant(0.0);
    //}

    jit_ins_branch( label_theta_exit );
    jit_ins_label( cont_1 );
    jit_label_t cont_2;
    jit_ins_branch( cont_2 , jit_ins_ge( eps.get_val() , jit_value( 1.0e-3 ) ) );
    // else if ( eps < 1.0e-3 ) {
    // ===============================================================================
    // Corner Case 3: c0->c0max even though c1 may be actually quite reasonable.
    // The ratio |c0|/c0max -> 1 but is still less than one, so that a 
    // series expansion is possible.
    // SERIES of acos(1-epsilon): Good to O(eps^6) or with this cutoff to O(10^{-18}) Computed with Maple.
    //  BTW: 1-epsilon = 1 - (c0max-c0abs)/c0max = 1-(1 - c0abs/c0max) = +c0abs/c0max
    //
    // ===============================================================================
    WordREG<REAL> sqtwo = sqrt( jit_constant(2.0) );
	      
    theta = 
      sqtwo * 
      sqrt(eps) * 
      ( jit_constant(1.0) + 
	( jit_constant(1/(REAL)12) + 
	  ( jit_constant(3/(REAL)160) + 
	    ( jit_constant(5/(REAL)896) + 
	      ( jit_constant(35/(REAL)18432) + 
		jit_constant(63/(REAL)90112) * eps ) * 
	      eps) *
	    eps) *
	  eps) *
	eps);

    jit_ins_branch( label_theta_exit );
    //} 
    jit_ins_label( cont_2 );
    //else {  
    // 
    theta = acos( c0abs/c0max );
    //}

    jit_ins_label( label_theta_exit );
      
      multi1d<WordREG<REAL> > f_site_re(3);
      multi1d<WordREG<REAL> > f_site_im(3);
	  
      multi1d<WordREG<REAL> > b1_site_re(3);
      multi1d<WordREG<REAL> > b1_site_im(3);
	      
      multi1d<WordREG<REAL> > b2_site_re(3);
      multi1d<WordREG<REAL> > b2_site_im(3);
		  
	  
	  
      WordREG<REAL> u = sqrt(c1/jit_constant(3.0))*cos(theta/jit_constant(3.0));
      WordREG<REAL> w = sqrt(c1)*sin(theta/jit_constant(3.0));
	  
      WordREG<REAL> u_sq = u*u;
      WordREG<REAL> w_sq = w*w;
	  
      WordREG<REAL> xi0,xi1;


      {
	jit_label_t label_90;
	jit_label_t cont_4;
	jit_value Nw_smallP = jit_ins_ge( (fabs( w )).get_val() , jit_value( 0.05 ) );
	jit_ins_branch( label_90 , Nw_smallP );
	{

	  xi0 = 
	    jit_constant(1.0) - 
	    (jit_constant(1.0/6.0)*w_sq*( jit_constant(1.0) - 
					  (jit_constant(1.0/20.)*w_sq*( jit_constant(1.0) - 
									(jit_constant(1.0/42.0)*w_sq ) ))));

	  jit_ins_branch( cont_4 );
	}
	jit_ins_label( label_90 );
	{
	  xi0 = sin(w)/w;
	}
	jit_ins_label( cont_4 );

	jit_label_t cont_3;
	jit_ins_branch( cont_3 , r_nobs );
	{
	  jit_label_t label_91;
	  jit_label_t cont_5;
	  jit_ins_branch( label_91 , Nw_smallP );
	  {
	    xi1 = 
	      jit_constant(-1.0)*
	      ( jit_constant((REAL)1/(REAL)3) - 
		jit_constant((REAL)1/(REAL)30)*w_sq*( jit_constant((REAL)1) - 
						      jit_constant((REAL)1/(REAL)28)*w_sq*( jit_constant((REAL)1) - 
											    jit_constant((REAL)1/(REAL)54)*w_sq ) ) );
	    jit_ins_branch( cont_5 );
	  }
	  {
	    jit_ins_label( label_91 );
	    xi1 = cos(w)/w_sq - sin(w)/(w_sq*w);
	  }
	  jit_ins_label( cont_5 );
	  jit_ins_label( cont_3 );
	}
      }      


      WordREG<REAL> cosu = cos(u);
      WordREG<REAL> sinu = sin(u);
      WordREG<REAL> cosw = cos(w);
      WordREG<REAL> sinw = sin(w);
      WordREG<REAL> sin2u = sin(jit_constant(2.0)*u);
      WordREG<REAL> cos2u = cos(jit_constant(2.0)*u);
      WordREG<REAL> ucosu = u*cosu;
      WordREG<REAL> usinu = u*sinu;
      WordREG<REAL> ucos2u = u*cos2u;
      WordREG<REAL> usin2u = u*sin2u;
	  
      WordREG<REAL> denum = jit_constant(9.0) * u_sq - w_sq;


      {
	WordREG<REAL> subexp1 = u_sq - w_sq;
	WordREG<REAL> subexp2 = jit_constant(8.0)*u_sq*cosw;
	WordREG<REAL> subexp3 = (jit_constant(3.0)*u_sq + w_sq)*xi0;
	    
	f_site_re[0] = ( (subexp1)*cos2u + cosu*subexp2 + jit_constant(2.0)*usinu*subexp3 ) / denum ;
	f_site_im[0] = ( (subexp1)*sin2u - sinu*subexp2 + jit_constant(2.0)*ucosu*subexp3 ) / denum ;
      }

      {
	WordREG<REAL> subexp = (jit_constant(3.0)*u_sq -w_sq)*xi0;
	    
	f_site_re[1] = (jit_constant(2.0)*(ucos2u - ucosu*cosw)+subexp*sinu)/denum;
	f_site_im[1] = (jit_constant(2.0)*(usin2u + usinu*cosw)+subexp*cosu)/denum;
      }

	  
      {
	WordREG<REAL> subexp=jit_constant(3.0)*xi0;
	    
	f_site_re[2] = (cos2u - cosu*cosw -usinu*subexp) /denum ;
	f_site_im[2] = (sin2u + sinu*cosw -ucosu*subexp) /denum ;
      }

      jit_label_t cont_6;
      jit_ins_branch( cont_6 , r_nobs );	  
      //if( dobs == true ) 
        {
	  multi1d<WordREG<REAL> > r_1_re(3);
	  multi1d<WordREG<REAL> > r_1_im(3);
	  multi1d<WordREG<REAL> > r_2_re(3);
	  multi1d<WordREG<REAL> > r_2_im(3);
	      
	  //	  r_1[0]=Double(2)*cmplx(u, u_sq-w_sq)*exp2iu
	  //          + 2.0*expmiu*( cmplx(8.0*u*cosw, -4.0*u_sq*cosw)
	  //	      + cmplx(u*(3.0*u_sq+w_sq),9.0*u_sq+w_sq)*xi0 );
	  {
	    WordREG<REAL> subexp1 = u_sq - w_sq;
	    WordREG<REAL> subexp2 =  jit_constant(8.0)*cosw + (jit_constant(3.0)*u_sq + w_sq)*xi0 ;
	    WordREG<REAL> subexp3 =  jit_constant(4.0)*u_sq*cosw - (jit_constant(9.0)*u_sq + w_sq)*xi0 ;
		
	    r_1_re[0] = jit_constant(2.0)*(ucos2u - sin2u *(subexp1)+ucosu*( subexp2 )- sinu*( subexp3 ) );
	    r_1_im[0] = jit_constant(2.0)*(usin2u + cos2u *(subexp1)-usinu*( subexp2 )- cosu*( subexp3 ) );
	  }
	      
	  // r_1[1]=cmplx(2.0, 4.0*u)*exp2iu + expmiu*cmplx(-2.0*cosw-(w_sq-3.0*u_sq)*xi0,2.0*u*cosw+6.0*u*xi0);
	  {
	    WordREG<REAL> subexp1 = cosw + jit_constant(3.0) * xi0;
	    WordREG<REAL> subexp2 = jit_constant(2.0)*cosw + xi0*(w_sq - jit_constant(3.0)*u_sq);
		
	    r_1_re[1] = jit_constant(2.0)*((cos2u - jit_constant(2.0)*usin2u) + usinu*( subexp1 )) - cosu*( subexp2 );
	    r_1_im[1] = jit_constant(2.0)*((sin2u + jit_constant(2.0)*ucos2u) + ucosu*( subexp1 )) + sinu*( subexp2 );
	  }
	      
	      
	  // r_1[2]=2.0*timesI(exp2iu)  +expmiu*cmplx(-3.0*u*xi0, cosw-3*xi0);
	  {
	    WordREG<REAL> subexp = cosw - jit_constant(3.0)*xi0;
	    r_1_re[2] = -jit_constant(2.0)*sin2u -jit_constant(3.0)*ucosu*xi0 + sinu*( subexp );
	    r_1_im[2] = jit_constant(2.0)*cos2u  +jit_constant(3.0)*usinu*xi0 + cosu*( subexp );
	  }
	      
	      
	  //r_2[0]=-2.0*exp2iu + 2*cmplx(0,u)*expmiu*cmplx(cosw+xi0+3*u_sq*xi1,
	  //						 4*u*xi0);
	  {
	    WordREG<REAL> subexp = cosw + xi0 + jit_constant(3.0)*u_sq*xi1;
	    r_2_re[0] = -jit_constant(2.0)*(cos2u + u*( jit_constant(4.0)*ucosu*xi0 - sinu*(subexp )) );
	    r_2_im[0] = -jit_constant(2.0)*(sin2u - u*( jit_constant(4.0)*usinu*xi0 + cosu*(subexp )) );
	  }
	      
	      
	  // r_2[1]= expmiu*cmplx(cosw+xi0-3.0*u_sq*xi1, 2.0*u*xi0);
	  // r_2[1] = timesMinusI(r_2[1]);
	  {
	    WordREG<REAL> subexp =  cosw + xi0 - jit_constant(3.0)*u_sq*xi1;
	    r_2_re[1] =  jit_constant(2.0)*ucosu*xi0 - sinu*( subexp ) ;
	    r_2_im[1] = jit_constant(-2.0)*usinu*xi0 - cosu*( subexp ) ;
	  }
	      
	  //r_2[2]=expmiu*cmplx(xi0, -3.0*u*xi1);
	  {
	    WordREG<REAL> subexp = jit_constant(3.0)*xi1;
		
	    r_2_re[2] =    cosu*xi0 - usinu*subexp ;
	    r_2_im[2] = -( sinu*xi0 + ucosu*subexp ) ;
	  }      
	      
	  WordREG<REAL> b_denum=jit_constant(2.0)*denum*denum;
	      
	      
	  for(int j=0; j < 3; j++) { 
		
	    {
	      WordREG<REAL> subexp1 = jit_constant(2.0)*u;
	      WordREG<REAL> subexp2 = jit_constant(3.0)*u_sq - w_sq;
	      WordREG<REAL> subexp3 = jit_constant(2.0)*(jit_constant(15.0)*u_sq + w_sq);
		  
	      b1_site_re[j]=( subexp1*r_1_re[j] + subexp2*r_2_re[j] - subexp3*f_site_re[j] )/b_denum;
	      b1_site_im[j]=( subexp1*r_1_im[j] + subexp2*r_2_im[j] - subexp3*f_site_im[j] )/b_denum;
	    }
		
	    { 
	      WordREG<REAL> subexp1 = jit_constant(3.0)*u;
	      WordREG<REAL> subexp2 = jit_constant(24.0)*u;

	      b2_site_re[j]=( r_1_re[j]- subexp1*r_2_re[j] - subexp2 * f_site_re[j] )/b_denum;
	      b2_site_im[j]=( r_1_im[j] -subexp1*r_2_im[j] - subexp2 * f_site_im[j] )/b_denum;
	    }
	  }

	  // Now flip the coefficients of the b-s
	  jit_label_t cont_7;
	  jit_ins_branch( cont_7 , jit_ins_not(c0_negativeP) );
	  //if( c0_negativeP ) 
	    {
	      //b1_site[0] = conj(b1_site[0]);
	      b1_site_im[0] *= jit_constant(-1.0);
		
	      //b1_site[1] = -conj(b1_site[1]);
	      b1_site_re[1] *= jit_constant(-1.0);
		
	      //b1_site[2] = conj(b1_site[2]);
	      b1_site_im[2] *= jit_constant(-1.0);
		
	      //b2_site[0] = -conj(b2_site[0]);
	      b2_site_re[0] *= jit_constant(-1.0);
		
	      //b2_site[1] = conj(b2_site[1]);
	      b2_site_im[1] *= jit_constant(-1.0);
		
	      //b2_site[2] = -conj(b2_site[2]);
	      b2_site_re[2] *= jit_constant(-1.0);
	    }
	  jit_ins_label( cont_7 );
	      
	  // Load back into the lattice sized object
	  //for(int j=0; j < 3; j++) {
		
	  b10_j.elem().elem().real() = b1_site_re[0];
	  b10_j.elem().elem().imag() = b1_site_im[0];
		
	  b20_j.elem().elem().real() = b2_site_re[0]; //ok
	  b20_j.elem().elem().imag() = b2_site_im[0]; //ok

	  b11_j.elem().elem().real() = b1_site_re[1]; //ok
	  b11_j.elem().elem().imag() = b1_site_im[1]; //ok
		
	  b21_j.elem().elem().real() = b2_site_re[1];
	  b21_j.elem().elem().imag() = b2_site_im[1];

	  b12_j.elem().elem().real() = b1_site_re[2];
	  b12_j.elem().elem().imag() = b1_site_im[2];
		
	  b22_j.elem().elem().real() = b2_site_re[2];
	  b22_j.elem().elem().imag() = b2_site_im[2];


	    //}
	} // end of if (dobs==true)
	jit_ins_label( cont_6 );

      // Now when everything is done flip signs of the b-s (can't do this before
      // as the unflipped f-s are needed to find the b-s
	  
	jit_label_t cont_8;
	jit_ins_branch( cont_8 , jit_ins_not(c0_negativeP) );
	//      if( c0_negativeP ) {
	    
	// f_site[0] = conj(f_site[0]);
	f_site_im[0] *= jit_constant(-1.0);
	    
	//f_site[1] = -conj(f_site[1]);
	f_site_re[1] *= jit_constant(-1.0);
	    
	//f_site[2] = conj(f_site[2]);
	f_site_im[2] *= jit_constant(-1.0);
	    
	//}
	jit_ins_label( cont_8 );
	
      // Load back into the lattice sized object
	//      for(int j=0; j < 3; j++) { 
	f0_j.elem().elem().real() = f_site_re[0];
	f0_j.elem().elem().imag() = f_site_im[0];

	f1_j.elem().elem().real() = f_site_re[1];
	f1_j.elem().elem().imag() = f_site_im[1];

	f2_j.elem().elem().real() = f_site_re[2];
	f2_j.elem().elem().imag() = f_site_im[2];
	//}
	jit_ins_label(label_exit);
    } // End of if( corner_caseP ) else {}

  return jit_get_cufunction("ptx_get_fs_bs.ptx");
}

#endif
