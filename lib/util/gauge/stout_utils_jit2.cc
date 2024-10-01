#include "qdp.h"

#include "chroma_config.h"

#if defined (QDP_IS_QDPJIT2)

using namespace QDP;


void function_get_fs_bs_exec(JitFunction& function, 
			     const LatticeColorMatrix& Q,
			     const LatticeColorMatrix& QQ,
			     multi1d<LatticeComplex>& f,
			     multi1d<LatticeComplex>& b1,
			     multi1d<LatticeComplex>& b2,
			     bool dobs)
{
  //QDPIO::cout << __FILE__ << ":" << __LINE__ << "\n";

  AddressLeaf addr_leaf(all);

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

  int th_count = Layout::sitesOnNode();
  
  WorkgroupGuardExec workgroupGuardExec(th_count , MG::get(Q.get_layout_ref()).sitesOnNode() );

  JitParam jit_dobs( QDP_get_global_cache().addJitParamBool( dobs ) );
  
  std::vector<QDPCache::ArgKey> ids;
  workgroupGuardExec.check(ids);
  ids.push_back( jit_dobs.get_id() );
  for(unsigned i=0; i < addr_leaf.ids.size(); ++i) 
    ids.push_back( addr_leaf.ids[i] );
  
  jit_launch(function,th_count,ids);
}







void function_get_fs_bs_build(JitFunction& function,
			      const LatticeColorMatrix& Q,
			      const LatticeColorMatrix& QQ,
			      multi1d<LatticeComplex>& f,
			      multi1d<LatticeComplex>& b1,
			      multi1d<LatticeComplex>& b2)
{
  llvm_start_new_function("get_fs_bs",__PRETTY_FUNCTION__);

  WorkgroupGuard workgroupGuard;

  ParamRef p_dobs   = llvm_add_param<bool>();

  ParamLeaf param_leaf(workgroupGuard);
  
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

  llvm::Value* r_idx = llvm_thread_idx();
  workgroupGuard.check(r_idx);

  llvm::Value*  r_dobs   = llvm_derefParam( p_dobs );
      

  auto Q_j  = Q_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto QQ_j = QQ_jit.elem(JitDeviceLayout::Coalesced,r_idx);

  auto f0_j = f0_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto f1_j = f1_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto f2_j = f2_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  
  auto b10_j = b10_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto b11_j = b11_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto b12_j = b12_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto b20_j = b20_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto b21_j = b21_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto b22_j = b22_jit.elem(JitDeviceLayout::Coalesced,r_idx);


  // Get the traces
  auto Q_site = Q_j.elem();   // col
  auto QQ_site = QQ_j.elem(); // col
  auto QQQ = QQ_site*Q_site;  // col
	

  auto trQQQ = realTrace(QQQ);  // col
  auto trQQ  = realTrace(QQ_site); // col

  auto c0    = 1./3. * trQQQ.elem().elem();  // eq 13
  auto c1    = 1./2. * trQQ.elem().elem();	 // eq 15 

  auto zero_word = create_word<REAL>(0.);
  auto m1_word = create_word<REAL>(-1.);
  auto onedot5 = create_word<REAL>(1.5);
  
  JitIf c1_lt_0p004( c1 < 4.0e-3 );  //    if( c1 < 4.0e-3  ) 
  {
    f0_j.elem().elem().real() = 1.0-c0*c0/720.0;
    f0_j.elem().elem().imag() =  -(c0/6.0)*(1.0-(c1/20.0)*(1.0-(c1/42.0))) ;
    
    f1_j.elem().elem().real() =  c0/24.0*(1.0-c1/15.0*(1.0-3.0*c1/112.0)) ;
    f1_j.elem().elem().imag() =  1.0-c1/6.0*(1.0-c1/20.0*(1.0-c1/42.0))-c0*c0/5040.0 ;
	  
    f2_j.elem().elem().real() = 0.5*(-1.0+c1/12.0*(1.0-c1/30.0*(1.0-c1/56.0))+c0*c0/20160.0);
    f2_j.elem().elem().imag() = 0.5*(c0/60.0*(1.0-c1/21.0*(1.0-c1/48.0)));

    
    JitIf if_dobs( r_dobs ); // dobs
    {
      //  partial f0/ partial c0
      b20_j.elem().elem().real() = -c0/360.0;
      b20_j.elem().elem().imag() =  -(1.0/6.0)*(1.0-(c1/20.0)*(1.0-c1/42.0));
	    
      // partial f0 / partial c1
      //
      b10_j.elem().elem().real() = zero_word;
      b10_j.elem().elem().imag() = (c0/120.0)*(1.0-c1/21.0);
	    
      // partial f1 / partial c0
      //
      b21_j.elem().elem().real() = (1.0/24.0)*(1.0-c1/15.0*(1.0-3.0*c1/112.0));
      b21_j.elem().elem().imag() = -c0/2520.0;
	    
	    
      // partial f1 / partial c1
      b11_j.elem().elem().real() = -c0/360.0*(1.0 - 3.0*c1/56.0 );
      b11_j.elem().elem().imag() = -1.0/6.0*(1.0-c1/10.0*(1.0-c1/28.0));
	    
      // partial f2/ partial c0
      b22_j.elem().elem().real() = 0.5*c0/10080.0;
      b22_j.elem().elem().imag() = 0.5*(  1.0/60.0*(1.0-c1/21.0*(1.0-c1/48.0)) );
	    
      // partial f2/ partial c1
      b12_j.elem().elem().real() = 0.5*(  1.0/12.0*(1.0-(2.0*c1/30.0)*(1.0-3.0*c1/112.0)) ); 
      b12_j.elem().elem().imag() = 0.5*( -c0/1260.0*(1.0-c1/24.0) );
    } // Dobs==true
    if_dobs.end();
  }
  c1_lt_0p004.els(); // if (c1 < 4.0e-3 )
  {
    auto c0abs = fabs( c0 );
    auto c0max = 2.0 * pow( c1 / 3.0 , onedot5 );
    auto eps = ( c0max - c0abs ) / c0max ;

    auto theta = create_word<REAL>();

    JitIf eps_lt_0( eps < 0.0 ); // epsilon < 0
    {
      theta = 0.0;
    }
    eps_lt_0.els();
    {
      JitIf eps_lt_0p001( eps < 0.001 ); // epsilon < 1.e-3
      {
	auto two = create_word<REAL>(2.0);
	auto sqtwo = sqrt( two );

	theta = sqtwo*sqrt(eps)*( 1.0 + ( (1/(REAL)12) + ( (3/(REAL)160) + ( (5/(REAL)896) + ( (35/(REAL)18432) + (63/(REAL)90112)*eps ) *eps) *eps) *eps) *eps);
      }
      eps_lt_0p001.els();                       // else
      {
	theta = acos( c0abs/c0max );
      }
      eps_lt_0p001.end();
    }
    eps_lt_0.end();
    

	

    auto f_site_re_0 = create_word<REAL>();
    auto f_site_re_1 = create_word<REAL>();
    auto f_site_re_2 = create_word<REAL>();
    
    auto f_site_im_0 = create_word<REAL>();
    auto f_site_im_1 = create_word<REAL>();
    auto f_site_im_2 = create_word<REAL>();

    ////
    
    auto b1_site_re_0 = create_word<REAL>();
    auto b1_site_re_1 = create_word<REAL>();
    auto b1_site_re_2 = create_word<REAL>();

    auto b1_site_im_0 = create_word<REAL>();
    auto b1_site_im_1 = create_word<REAL>();
    auto b1_site_im_2 = create_word<REAL>();
    
    auto b2_site_re_0 = create_word<REAL>();
    auto b2_site_re_1 = create_word<REAL>();
    auto b2_site_re_2 = create_word<REAL>();

    auto b2_site_im_0 = create_word<REAL>();
    auto b2_site_im_1 = create_word<REAL>();
    auto b2_site_im_2 = create_word<REAL>();

    ////
	  
    auto u = sqrt(c1/3)*cos(theta/3);
    auto w = sqrt(c1)*sin(theta/3);
	  
    auto u_sq = u*u;
    auto w_sq = w*w;

    
    auto xi0 = create_word< REAL >();
    auto xi1 = create_word< REAL >(0.0);

    JitIf w_small( ( fabs( w ) < 0.05 ) );
    {
      xi0 = (REAL)1 - ((REAL)1/(REAL)6)*w_sq*( 1 - ((REAL)1/(REAL)20)*w_sq*( (REAL)1 - ((REAL)1/(REAL)42)*w_sq ) );
    }
    w_small.els();
    {
      xi0 = sin(w)/w;
    }
    w_small.end();

    
    JitIf if_dobs2( r_dobs );
    {
      JitIf w_small( fabs( w ) < 0.05 );
      {
	xi1 = -1*( ((REAL)1/(REAL)3) - ((REAL)1/(REAL)30)*w_sq*( (REAL)1 - ((REAL)1/(REAL)28)*w_sq*( (REAL)1 - ((REAL)1/(REAL)54)*w_sq ) ) );
      }
      w_small.els();
      {
	xi1 = cos(w)/w_sq - sin(w)/(w_sq*w);
      }
      w_small.end();
    }
    if_dobs2.end();



    auto cosu = cos(u);
    auto sinu = sin(u);
    auto cosw = cos(w);
    auto sinw = sin(w);
    auto sin2u = sin(2.0*u);
    auto cos2u = cos(2.0*u);
    auto ucosu = u*cosu;
    auto usinu = u*sinu;
    auto ucos2u = u*cos2u;
    auto usin2u = u*sin2u;
	  
    auto denum = 9.0 * u_sq - w_sq;

    {
      auto subexp1 = u_sq - w_sq;
      auto subexp2 = 8.0*u_sq*cosw;
      auto subexp3 = (3.0*u_sq + w_sq)*xi0;
	    
      f_site_re_0 = ( (subexp1)*cos2u + cosu*subexp2 + 2.0*usinu*subexp3 ) / denum ;
      f_site_im_0 = ( (subexp1)*sin2u - sinu*subexp2 + 2.0*ucosu*subexp3 ) / denum ;
    }

    {
      auto subexp = (3.0*u_sq -w_sq)*xi0;
	    
      f_site_re_1 = (2.0*(ucos2u - ucosu*cosw)+subexp*sinu)/denum;
      f_site_im_1 = (2.0*(usin2u + usinu*cosw)+subexp*cosu)/denum;
    }

	  
    {
      auto subexp=3.0*xi0;
	    
      f_site_re_2 = (cos2u - cosu*cosw -usinu*subexp) /denum ;
      f_site_im_2 = (sin2u + sinu*cosw -ucosu*subexp) /denum ;
    }

    
    JitIf if_dobs3( r_dobs );
    {
      auto r_1_re_0 = create_word<REAL>();
      auto r_1_re_1 = create_word<REAL>();
      auto r_1_re_2 = create_word<REAL>();

      auto r_1_im_0 = create_word<REAL>();
      auto r_1_im_1 = create_word<REAL>();
      auto r_1_im_2 = create_word<REAL>();

      auto r_2_re_0 = create_word<REAL>();
      auto r_2_re_1 = create_word<REAL>();
      auto r_2_re_2 = create_word<REAL>();

      auto r_2_im_0 = create_word<REAL>();
      auto r_2_im_1 = create_word<REAL>();
      auto r_2_im_2 = create_word<REAL>();

      //	  r_1[0]=Double(2)*cmplx(u, u_sq-w_sq)*exp2iu
      //          + 2.0*expmiu*( cmplx(8.0*u*cosw, -4.0*u_sq*cosw)
      //	      + cmplx(u*(3.0*u_sq+w_sq),9.0*u_sq+w_sq)*xi0 );
      {
	auto subexp1 = u_sq - w_sq;
	auto subexp2 = 8.0*cosw + (3.0*u_sq + w_sq)*xi0 ;
	auto subexp3 = 4.0*u_sq*cosw - (9.0*u_sq + w_sq)*xi0 ;
		
	r_1_re_0 = 2.0*(ucos2u - sin2u *(subexp1)+ucosu*( subexp2 )- sinu*( subexp3 ) );
	r_1_im_0 = 2.0*(usin2u + cos2u *(subexp1)-usinu*( subexp2 )- cosu*( subexp3 ) );
      }
	      
      // r_1[1]=cmplx(2.0, 4.0*u)*exp2iu + expmiu*cmplx(-2.0*cosw-(w_sq-3.0*u_sq)*xi0,2.0*u*cosw+6.0*u*xi0);
      {
	auto subexp1 = cosw + 3.0 * xi0;
	auto subexp2 = 2.0*cosw + xi0*(w_sq - 3.0*u_sq);
		
	r_1_re_1 = 2.0*((cos2u - 2.0*usin2u) + usinu*( subexp1 )) - cosu*( subexp2 );
	r_1_im_1 = 2.0*((sin2u + 2.0*ucos2u) + ucosu*( subexp1 )) + sinu*( subexp2 );
      }
	      
	      
      // r_1[2]=2.0*timesI(exp2iu)  +expmiu*cmplx(-3.0*u*xi0, cosw-3*xi0);
      {
	auto subexp = cosw - 3.0*xi0;
	r_1_re_2 = -2.0*sin2u -3.0*ucosu*xi0 + sinu*( subexp );
	r_1_im_2 = 2.0*cos2u  +3.0*usinu*xi0 + cosu*( subexp );
      }
	      
	      
      //r_2[0]=-2.0*exp2iu + 2*cmplx(0,u)*expmiu*cmplx(cosw+xi0+3*u_sq*xi1,
      //						 4*u*xi0);
      {
	auto subexp = cosw + xi0 + 3.0*u_sq*xi1;
	r_2_re_0 = -2.0*(cos2u + u*( 4.0*ucosu*xi0 - sinu*(subexp )) );
	r_2_im_0 = -2.0*(sin2u - u*( 4.0*usinu*xi0 + cosu*(subexp )) );
      }
	      
	      
      // r_2[1]= expmiu*cmplx(cosw+xi0-3.0*u_sq*xi1, 2.0*u*xi0);
      // r_2[1] = timesMinusI(r_2[1]);
      {
	auto subexp =  cosw + xi0 - 3.0*u_sq*xi1;
	r_2_re_1 = 2.0*ucosu*xi0 - sinu*( subexp ) ;
	r_2_im_1 = -2.0*usinu*xi0 - cosu*( subexp ) ;
      }
	      
      //r_2[2]=expmiu*cmplx(xi0, -3.0*u*xi1);
      {
	auto subexp = 3.0*xi1;
		
	r_2_re_2 =    cosu*xi0 - usinu*subexp ;
	r_2_im_2 = -( sinu*xi0 + ucosu*subexp ) ;
      }      
	      
      auto b_denum=2.0*denum*denum;
	      
      auto subexp1_b1 = 2.0*u;
      auto subexp2_b1 = 3.0*u_sq - w_sq;
      auto subexp3_b1 = 2.0*(15.0*u_sq + w_sq);
      
      auto subexp1_b2 = 3.0*u;
      auto subexp2_b2 = 24.0*u;

      b1_site_re_0 = ( subexp1_b1*r_1_re_0 + subexp2_b1*r_2_re_0 - subexp3_b1*f_site_re_0 )/b_denum;
      b1_site_re_1 = ( subexp1_b1*r_1_re_1 + subexp2_b1*r_2_re_1 - subexp3_b1*f_site_re_1 )/b_denum;
      b1_site_re_2 = ( subexp1_b1*r_1_re_2 + subexp2_b1*r_2_re_2 - subexp3_b1*f_site_re_2 )/b_denum;
      
      b1_site_im_0 = ( subexp1_b1*r_1_im_0 + subexp2_b1*r_2_im_0 - subexp3_b1*f_site_im_0 )/b_denum;
      b1_site_im_1 = ( subexp1_b1*r_1_im_1 + subexp2_b1*r_2_im_1 - subexp3_b1*f_site_im_1 )/b_denum;
      b1_site_im_2 = ( subexp1_b1*r_1_im_2 + subexp2_b1*r_2_im_2 - subexp3_b1*f_site_im_2 )/b_denum;

      b2_site_re_0=( r_1_re_0- subexp1_b2*r_2_re_0 - subexp2_b2 * f_site_re_0 )/b_denum;
      b2_site_re_1=( r_1_re_1- subexp1_b2*r_2_re_1 - subexp2_b2 * f_site_re_1 )/b_denum;
      b2_site_re_2=( r_1_re_2- subexp1_b2*r_2_re_2 - subexp2_b2 * f_site_re_2 )/b_denum;
      
      b2_site_im_0=( r_1_im_0 -subexp1_b2*r_2_im_0 - subexp2_b2 * f_site_im_0 )/b_denum;
      b2_site_im_1=( r_1_im_1 -subexp1_b2*r_2_im_1 - subexp2_b2 * f_site_im_1 )/b_denum;
      b2_site_im_2=( r_1_im_2 -subexp1_b2*r_2_im_2 - subexp2_b2 * f_site_im_2 )/b_denum;

      
      // Now flip the coefficients of the b-s

      JitIf if_c0_neg( c0 < 0.0 );   // if( c0_negativeP ) 
      {
	b1_site_im_0 *= m1_word;
	b1_site_re_1 *= m1_word;
	b1_site_im_2 *= m1_word;
	b2_site_re_0 *= m1_word;
	b2_site_im_1 *= m1_word;
	b2_site_re_2 *= m1_word;
      }
      if_c0_neg.end();

      b10_j.elem().elem().real() = b1_site_re_0;
      b10_j.elem().elem().imag() = b1_site_im_0;
		
      b20_j.elem().elem().real() = b2_site_re_0;
      b20_j.elem().elem().imag() = b2_site_im_0;

      b11_j.elem().elem().real() = b1_site_re_1;
      b11_j.elem().elem().imag() = b1_site_im_1;
		
      b21_j.elem().elem().real() = b2_site_re_1;
      b21_j.elem().elem().imag() = b2_site_im_1;

      b12_j.elem().elem().real() = b1_site_re_2;
      b12_j.elem().elem().imag() = b1_site_im_2;
		
      b22_j.elem().elem().real() = b2_site_re_2;
      b22_j.elem().elem().imag() = b2_site_im_2;
    } 
    if_dobs3.end();


    f0_j.elem().elem().real() = f_site_re_0;
    f0_j.elem().elem().imag() = f_site_im_0;

    f1_j.elem().elem().real() = f_site_re_1;
    f1_j.elem().elem().imag() = f_site_im_1;

    f2_j.elem().elem().real() = f_site_re_2;
    f2_j.elem().elem().imag() = f_site_im_2;

    JitIf if_c0_neg2( c0 < 0.0 );   // if( c0_negativeP ) 
    {
      f0_j.elem().elem().imag() *= m1_word;
      f1_j.elem().elem().real() *= m1_word;
      f2_j.elem().elem().imag() *= m1_word;
    }
    if_c0_neg2.end();

  }
  c1_lt_0p004.end(); // if (c1 < 4.0e-3 )

  jit_get_function(function);
}







#endif
