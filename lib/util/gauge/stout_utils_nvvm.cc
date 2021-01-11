#include "qdp.h"

#ifdef QDP_IS_QDPJIT
#ifdef QDPJIT_IS_QDPJITNVVM

using namespace QDP;

//#define QDP_JIT_NVVM_USE_LEGACY_LAUNCH

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

  int lo = 0;
  int hi = Layout::sitesOnNode();

  JitParam jit_lo( QDP_get_global_cache().addJitParamInt( lo ) );
  JitParam jit_hi( QDP_get_global_cache().addJitParamInt( hi ) );
  JitParam jit_dobs( QDP_get_global_cache().addJitParamBool( dobs ) );
  
  std::vector<QDPCache::ArgKey> ids;
  
  ids.push_back( jit_lo.get_id() );
  ids.push_back( jit_hi.get_id() );
  ids.push_back( jit_dobs.get_id() );
  
  for(unsigned i=0; i < addr_leaf.ids.size(); ++i) 
    ids.push_back( addr_leaf.ids[i] );
  
  jit_launch(function,hi-lo,ids);
}



using real_t = RScalarREG< WordREG<REAL> >;

namespace
{
  template <class T>
  decltype(auto) poke_real( T& t )
  {
    typedef typename JITType< real_t >::Type_t T_jit;
    return T_jit( t.elem().elem().real() );
  }

  template <class T>
  decltype(auto) poke_imag( T& t )
  {
    typedef typename JITType< real_t >::Type_t T_jit;
    return (T_jit( t.elem().elem().imag() ));
  }
}


void function_get_fs_bs_build(JitFunction& function,
			      const LatticeColorMatrix& Q,
			      const LatticeColorMatrix& QQ,
			      multi1d<LatticeComplex>& f,
			      multi1d<LatticeComplex>& b1,
			      multi1d<LatticeComplex>& b2)
{
  if (ptx_db::db_enabled)
    {
      llvm_ptx_db( function , __PRETTY_FUNCTION__ );
      if (!function.empty())
	return;
    }

  //std::cout << __PRETTY_FUNCTION__ << ": entering\n";


  llvm_start_new_function("get_fs_bs",__PRETTY_FUNCTION__);

  ParamRef  p_lo     = llvm_add_param<int>();
  ParamRef  p_hi     = llvm_add_param<int>();
  ParamRef  p_dobs   = llvm_add_param<bool>();

  ParamLeaf param_leaf;
  
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

  llvm::Value*  r_lo     = llvm_derefParam( p_lo );
  llvm::Value*  r_hi     = llvm_derefParam( p_hi );
  llvm::Value*  r_dobs   = llvm_derefParam( p_dobs );
  //llvm::Value*  r_nobs   = llvm_not( r_dobs );
      
  llvm::Value*  r_idx = llvm_thread_idx();
      
  llvm_cond_exit( llvm_ge( r_idx , r_hi ) );


  auto& Q_j  = Q_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto& QQ_j = QQ_jit.elem(JitDeviceLayout::Coalesced,r_idx);

  auto& f0_j = f0_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto& f1_j = f1_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto& f2_j = f2_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  
  auto& b10_j = b10_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto& b11_j = b11_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto& b12_j = b12_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto& b20_j = b20_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto& b21_j = b21_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto& b22_j = b22_jit.elem(JitDeviceLayout::Coalesced,r_idx);


  // Get the traces
  PColorMatrixREG< RComplexREG< WordREG<REAL> >, Nc>  Q_site = Q_j.elem();
  PColorMatrixREG< RComplexREG< WordREG<REAL> >, Nc>  QQ_site = QQ_j.elem();
  PColorMatrixREG< RComplexREG< WordREG<REAL> >, Nc>  QQQ = QQ_site*Q_site;
	

  PScalarREG< RScalarREG< WordREG<REAL> > > trQQQ = realTrace(QQQ);
  PScalarREG< RScalarREG< WordREG<REAL> > > trQQ  = realTrace(QQ_site);

  real_t c0    = real_t( 1./3.) * trQQQ.elem();  // eq 13
  real_t c1    = real_t( 1./2.) * trQQ.elem();	 // eq 15 


  auto lv_c1_lt_0p004 = real_t( c1 < real_t( 4.0e-3 ) ).elem().get_val();

  JitIf c1_lt_0p004( lv_c1_lt_0p004 );   //    if( c1 < 4.0e-3  ) 
  {

    poke_real(f0_j) = real_t(1.0) - c0 * c0 / real_t(720.0);
    poke_imag(f0_j) =  -( c0 / real_t(6.0) )*( real_t(1.0) -(c1/real_t(20.0))*(real_t(1.0)-(c1/real_t(42.0)))) ;

    poke_real(f1_j) =  c0/real_t(24.0)*(real_t(1.0)-c1/real_t(15.0)*(real_t(1.0)-real_t(3.0)*c1/real_t(112.0))) ;
    poke_imag(f1_j) =  real_t(1.0)-c1/real_t(6.0)*(real_t(1.0)-c1/real_t(20.0)*(real_t(1.0)-c1/real_t(42.0)))-c0*c0/real_t(5040.0);
	  
    poke_real(f2_j) = real_t(0.5)*(real_t(-1.0)+c1/real_t(12.0)*(real_t(1.0)-c1/real_t(30.0)*(real_t(1.0)-c1/real_t(56.0)))+c0*c0/real_t(20160.0));
    poke_imag(f2_j) = real_t(0.5)*(c0/real_t(60.0)*(real_t(1.0)-c1/real_t(21.0)*(real_t(1.0)-c1/real_t(48.0))));


    JitIf if_dobs( r_dobs ); // dobs
    {
      poke_real(b20_j) = -c0/real_t(360.0);
      poke_imag(b20_j) =  -real_t(1.0/6.0)*(real_t(1.0)-(c1/real_t(20.0))*(real_t(1.0)-c1/real_t(42.0)));
	    
      // partial f0 / partial c1
      //
      poke_real(b10_j) = real_t(0);
      poke_imag(b10_j) = (c0/real_t(120.0))*(real_t(1.0)-c1/real_t(21.0));
	    
      // partial f1 / partial c0
      //
      poke_real(b21_j) = real_t(1.0/24.0)*(real_t(1.0)-c1/real_t(15.0)*(real_t(1.0)-real_t(3.0)*c1/real_t(112.0)));
      poke_imag(b21_j) = -c0/real_t(2520.0);
	    
	    
      // partial f1 / partial c1
      poke_real(b11_j) = -c0/real_t(360.0)*(real_t(1.0) - real_t(3.0)*c1/real_t(56.0) );
      poke_imag(b11_j) = -real_t(1.0/6.0)*(real_t(1.0)-c1/real_t(10.0)*(real_t(1.0)-c1/real_t(28.0)));
	    
      // partial f2/ partial c0
      poke_real(b22_j) = real_t(0.5)*c0/real_t(10080.0);
      poke_imag(b22_j) = real_t(0.5)*(  real_t(1.0/60.0)*(real_t(1.0)-c1/real_t(21.0)*(real_t(1.0)-c1/real_t(48.0))) );
	    
      // partial f2/ partial c1
      poke_real(b12_j) = real_t(0.5)*(  real_t(1.0/12.0)*(real_t(1.0)-(real_t(2.0)*c1/real_t(30.0))*(real_t(1.0)-real_t(3.0)*c1/real_t(112.0))) ); 
      poke_imag(b12_j) = real_t(0.5)*( -c0/real_t(1260.0)*(real_t(1.0)-c1/real_t(24.0)) );

    } // Dobs==true
    if_dobs.end();
    
  }
  c1_lt_0p004.els(); // if (c1 < 4.0e-3 )
  {
    
    real_t c0abs = fabs( c0 );
    real_t c0max = real_t( 2.0 ) * pow( c1 / real_t( 3.0 ) , real_t( 1.5 ) );
    real_t theta;
    real_t eps = ( c0max - c0abs ) / c0max ;


    auto theta_stack = stack_alloc< real_t >();
    
    auto lv_eps_lt_0 = real_t(eps < real_t( 0.0 )).elem().get_val();
    
    JitIf eps_lt_0( lv_eps_lt_0 ); // epsilon < 0
    {
      theta_stack = real_t( 0.0 );
    }
    eps_lt_0.els();
    {
      auto lv_eps_lt_0p001 = real_t(eps < real_t( 0.001 )).elem().get_val();

      JitIf eps_lt_0p001( lv_eps_lt_0p001 ); // epsilon < 1.e-3
      {
	real_t sqtwo = sqrt( real_t(2.0) );
	real_t theta_tmp;
	theta_tmp = 
	  sqtwo * 
	  sqrt(eps) * 
	  ( real_t(1.0) + 
	    ( real_t(1/(REAL)12) + 
	      ( real_t(3/(REAL)160) + 
		( real_t(5/(REAL)896) + 
		  ( real_t(35/(REAL)18432) + 
		    real_t(63/(REAL)90112) * eps ) * 
		  eps) *
		eps) *
	      eps) *
	    eps);
	theta_stack = theta_tmp;
      }
      eps_lt_0p001.els();                       // else
      {
	real_t theta_tmp2;
	theta_tmp2 = acos( c0abs/c0max );
	theta_stack = theta_tmp2;
      }
      eps_lt_0p001.end();
    }
    eps_lt_0.end();
    
    theta.setup( theta_stack );


	
    multi1d< real_t > f_site_re(3);
    multi1d< real_t > f_site_im(3);
	  
    auto b1_site_re_0_stack = stack_alloc< real_t >();
    auto b1_site_re_1_stack = stack_alloc< real_t >();
    auto b1_site_re_2_stack = stack_alloc< real_t >();

    auto b1_site_im_0_stack = stack_alloc< real_t >();
    auto b1_site_im_1_stack = stack_alloc< real_t >();
    auto b1_site_im_2_stack = stack_alloc< real_t >();
    
    auto b2_site_re_0_stack = stack_alloc< real_t >();
    auto b2_site_re_1_stack = stack_alloc< real_t >();
    auto b2_site_re_2_stack = stack_alloc< real_t >();

    auto b2_site_im_0_stack = stack_alloc< real_t >();
    auto b2_site_im_1_stack = stack_alloc< real_t >();
    auto b2_site_im_2_stack = stack_alloc< real_t >();
	  
	  
    real_t u = sqrt(c1/real_t(3.0))*cos(theta/real_t(3.0));
    real_t w = sqrt(c1)*sin(theta/real_t(3.0));
	  
    real_t u_sq = u*u;
    real_t w_sq = w*w;
	  
    real_t xi0,xi1;

    { // xi0

      auto xi0_stack = stack_alloc< real_t >();

      llvm::Value* lv_w_small = real_t( fabs( w ) < real_t( 0.05 ) ).elem().get_val();
	
      JitIf w_small( lv_w_small );
      {
	real_t xi0_tmp0 =
	  real_t(1.0) - 
	  (real_t(1.0/6.0)*w_sq*( real_t(1.0) - 
				  (real_t(1.0/20.)*w_sq*( real_t(1.0) - 
							  (real_t(1.0/42.0)*w_sq ) ))));
	xi0_stack = xi0_tmp0;
      }
      w_small.els();
      {
	real_t xi0_tmp1 = sin(w)/w;
	xi0_stack = xi0_tmp1;
      }
      w_small.end();
      
      xi0.setup( xi0_stack );
    } // xi0

    
    auto xi1_stack = stack_alloc< real_t >();

    xi1_stack = real_t(0.0);
    
    JitIf if_dobs2( r_dobs );
    {
      { // xi1
	
	llvm::Value* lv_w_small = real_t( fabs( w ) < real_t( 0.05 ) ).elem().get_val();
	
	JitIf w_small( lv_w_small );
	{
	  real_t xi1_tmp0 = 	    
	    real_t(-1.0)*
	    ( real_t((REAL)1/(REAL)3) - 
	      real_t((REAL)1/(REAL)30)*w_sq*( real_t((REAL)1) - 
					      real_t((REAL)1/(REAL)28)*w_sq*( real_t((REAL)1) - 
									      real_t((REAL)1/(REAL)54)*w_sq ) ) );
	  xi1_stack = xi1_tmp0;
	}
	w_small.els();
	{
	  real_t xi1_tmp1 = cos(w)/w_sq - sin(w)/(w_sq*w);
	  xi1_stack = xi1_tmp1;
	}
	w_small.end();
      } // xi1
    }
    if_dobs2.end();

    xi1.setup( xi1_stack );

    real_t cosu = cos(u);
    real_t sinu = sin(u);
    real_t cosw = cos(w);
    real_t sinw = sin(w);
    real_t sin2u = sin(real_t(2.0)*u);
    real_t cos2u = cos(real_t(2.0)*u);
    real_t ucosu = u*cosu;
    real_t usinu = u*sinu;
    real_t ucos2u = u*cos2u;
    real_t usin2u = u*sin2u;
	  
    real_t denum = real_t(9.0) * u_sq - w_sq;

    {
      real_t subexp1 = u_sq - w_sq;
      real_t subexp2 = real_t(8.0)*u_sq*cosw;
      real_t subexp3 = (real_t(3.0)*u_sq + w_sq)*xi0;
	    
      f_site_re[0] = ( (subexp1)*cos2u + cosu*subexp2 + real_t(2.0)*usinu*subexp3 ) / denum ;
      f_site_im[0] = ( (subexp1)*sin2u - sinu*subexp2 + real_t(2.0)*ucosu*subexp3 ) / denum ;
    }

    {
      real_t subexp = (real_t(3.0)*u_sq -w_sq)*xi0;
	    
      f_site_re[1] = (real_t(2.0)*(ucos2u - ucosu*cosw)+subexp*sinu)/denum;
      f_site_im[1] = (real_t(2.0)*(usin2u + usinu*cosw)+subexp*cosu)/denum;
    }

	  
    {
      real_t subexp=real_t(3.0)*xi0;
	    
      f_site_re[2] = (cos2u - cosu*cosw -usinu*subexp) /denum ;
      f_site_im[2] = (sin2u + sinu*cosw -ucosu*subexp) /denum ;
    }

    
    JitIf if_dobs3( r_dobs );
    {
      multi1d<real_t > r_1_re(3);
      multi1d<real_t > r_1_im(3);
      multi1d<real_t > r_2_re(3);
      multi1d<real_t > r_2_im(3);
	      
      //	  r_1[0]=Double(2)*cmplx(u, u_sq-w_sq)*exp2iu
      //          + 2.0*expmiu*( cmplx(8.0*u*cosw, -4.0*u_sq*cosw)
      //	      + cmplx(u*(3.0*u_sq+w_sq),9.0*u_sq+w_sq)*xi0 );
      {
	real_t subexp1 = u_sq - w_sq;
	real_t subexp2 =  real_t(8.0)*cosw + (real_t(3.0)*u_sq + w_sq)*xi0 ;
	real_t subexp3 =  real_t(4.0)*u_sq*cosw - (real_t(9.0)*u_sq + w_sq)*xi0 ;
		
	r_1_re[0] = real_t(2.0)*(ucos2u - sin2u *(subexp1)+ucosu*( subexp2 )- sinu*( subexp3 ) );
	r_1_im[0] = real_t(2.0)*(usin2u + cos2u *(subexp1)-usinu*( subexp2 )- cosu*( subexp3 ) );
      }
	      
      // r_1[1]=cmplx(2.0, 4.0*u)*exp2iu + expmiu*cmplx(-2.0*cosw-(w_sq-3.0*u_sq)*xi0,2.0*u*cosw+6.0*u*xi0);
      {
	real_t subexp1 = cosw + real_t(3.0) * xi0;
	real_t subexp2 = real_t(2.0)*cosw + xi0*(w_sq - real_t(3.0)*u_sq);
		
	r_1_re[1] = real_t(2.0)*((cos2u - real_t(2.0)*usin2u) + usinu*( subexp1 )) - cosu*( subexp2 );
	r_1_im[1] = real_t(2.0)*((sin2u + real_t(2.0)*ucos2u) + ucosu*( subexp1 )) + sinu*( subexp2 );
      }
	      
	      
      // r_1[2]=2.0*timesI(exp2iu)  +expmiu*cmplx(-3.0*u*xi0, cosw-3*xi0);
      {
	real_t subexp = cosw - real_t(3.0)*xi0;
	r_1_re[2] = -real_t(2.0)*sin2u -real_t(3.0)*ucosu*xi0 + sinu*( subexp );
	r_1_im[2] = real_t(2.0)*cos2u  +real_t(3.0)*usinu*xi0 + cosu*( subexp );
      }
	      
	      
      //r_2[0]=-2.0*exp2iu + 2*cmplx(0,u)*expmiu*cmplx(cosw+xi0+3*u_sq*xi1,
      //						 4*u*xi0);
      {
	real_t subexp = cosw + xi0 + real_t(3.0)*u_sq*xi1;
	r_2_re[0] = -real_t(2.0)*(cos2u + u*( real_t(4.0)*ucosu*xi0 - sinu*(subexp )) );
	r_2_im[0] = -real_t(2.0)*(sin2u - u*( real_t(4.0)*usinu*xi0 + cosu*(subexp )) );
      }
	      
	      
      // r_2[1]= expmiu*cmplx(cosw+xi0-3.0*u_sq*xi1, 2.0*u*xi0);
      // r_2[1] = timesMinusI(r_2[1]);
      {
	real_t subexp =  cosw + xi0 - real_t(3.0)*u_sq*xi1;
	r_2_re[1] =  real_t(2.0)*ucosu*xi0 - sinu*( subexp ) ;
	r_2_im[1] = real_t(-2.0)*usinu*xi0 - cosu*( subexp ) ;
      }
	      
      //r_2[2]=expmiu*cmplx(xi0, -3.0*u*xi1);
      {
	real_t subexp = real_t(3.0)*xi1;
		
	r_2_re[2] =    cosu*xi0 - usinu*subexp ;
	r_2_im[2] = -( sinu*xi0 + ucosu*subexp ) ;
      }      
	      
      real_t b_denum=real_t(2.0)*denum*denum;
	      
      real_t subexp1_b1 = real_t(2.0)*u;
      real_t subexp2_b1 = real_t(3.0)*u_sq - w_sq;
      real_t subexp3_b1 = real_t(2.0)*(real_t(15.0)*u_sq + w_sq);
      
      real_t subexp1_b2 = real_t(3.0)*u;
      real_t subexp2_b2 = real_t(24.0)*u;

      b1_site_re_0_stack = ( subexp1_b1*r_1_re[0] + subexp2_b1*r_2_re[0] - subexp3_b1*f_site_re[0] )/b_denum;
      b1_site_re_1_stack = ( subexp1_b1*r_1_re[1] + subexp2_b1*r_2_re[1] - subexp3_b1*f_site_re[1] )/b_denum;
      b1_site_re_2_stack = ( subexp1_b1*r_1_re[2] + subexp2_b1*r_2_re[2] - subexp3_b1*f_site_re[2] )/b_denum;
      
      b1_site_im_0_stack = ( subexp1_b1*r_1_im[0] + subexp2_b1*r_2_im[0] - subexp3_b1*f_site_im[0] )/b_denum;
      b1_site_im_1_stack = ( subexp1_b1*r_1_im[1] + subexp2_b1*r_2_im[1] - subexp3_b1*f_site_im[1] )/b_denum;
      b1_site_im_2_stack = ( subexp1_b1*r_1_im[2] + subexp2_b1*r_2_im[2] - subexp3_b1*f_site_im[2] )/b_denum;

      b2_site_re_0_stack=( r_1_re[0]- subexp1_b2*r_2_re[0] - subexp2_b2 * f_site_re[0] )/b_denum;
      b2_site_re_1_stack=( r_1_re[1]- subexp1_b2*r_2_re[1] - subexp2_b2 * f_site_re[1] )/b_denum;
      b2_site_re_2_stack=( r_1_re[2]- subexp1_b2*r_2_re[2] - subexp2_b2 * f_site_re[2] )/b_denum;
      
      b2_site_im_0_stack=( r_1_im[0] -subexp1_b2*r_2_im[0] - subexp2_b2 * f_site_im[0] )/b_denum;
      b2_site_im_1_stack=( r_1_im[1] -subexp1_b2*r_2_im[1] - subexp2_b2 * f_site_im[1] )/b_denum;
      b2_site_im_2_stack=( r_1_im[2] -subexp1_b2*r_2_im[2] - subexp2_b2 * f_site_im[2] )/b_denum;

      
      // Now flip the coefficients of the b-s

      multi1d<real_t > b1_site_re_phi1(3);
      multi1d<real_t > b1_site_im_phi1(3);
	      
      multi1d<real_t > b2_site_re_phi1(3);
      multi1d<real_t > b2_site_im_phi1(3);

      auto lv_c0_neg = real_t( c0 < real_t(0.0) ).elem().get_val();
      
      JitIf if_c0_neg( lv_c0_neg );   // if( c0_negativeP ) 
      {
	b1_site_im_0_stack *= real_t(-1.0);
	b1_site_re_1_stack *= real_t(-1.0);
	b1_site_im_2_stack *= real_t(-1.0);
	b2_site_re_0_stack *= real_t(-1.0);
	b2_site_im_1_stack *= real_t(-1.0);
	b2_site_re_2_stack *= real_t(-1.0);
      }
      if_c0_neg.end();
  
      multi1d<real_t > f_site_re(3);
      multi1d<real_t > f_site_im(3);
	  
      multi1d<real_t > b1_site_re(3);
      multi1d<real_t > b1_site_im(3);
      
      multi1d<real_t > b2_site_re(3);
      multi1d<real_t > b2_site_im(3);

      b1_site_re[0].setup( b1_site_re_0_stack );
      b1_site_re[1].setup( b1_site_re_1_stack );
      b1_site_re[2].setup( b1_site_re_2_stack );

      b1_site_im[0].setup( b1_site_im_0_stack );
      b1_site_im[1].setup( b1_site_im_1_stack );
      b1_site_im[2].setup( b1_site_im_2_stack );

      b2_site_re[0].setup( b2_site_re_0_stack );
      b2_site_re[1].setup( b2_site_re_1_stack );
      b2_site_re[2].setup( b2_site_re_2_stack );

      b2_site_im[0].setup( b2_site_im_0_stack );
      b2_site_im[1].setup( b2_site_im_1_stack );
      b2_site_im[2].setup( b2_site_im_2_stack );

		
      poke_real(b10_j) = b1_site_re[0];
      poke_imag(b10_j) = b1_site_im[0];
		
      poke_real(b20_j) = b2_site_re[0];
      poke_imag(b20_j) = b2_site_im[0];

      poke_real(b11_j) = b1_site_re[1];
      poke_imag(b11_j) = b1_site_im[1];
		
      poke_real(b21_j) = b2_site_re[1];
      poke_imag(b21_j) = b2_site_im[1];

      poke_real(b12_j) = b1_site_re[2];
      poke_imag(b12_j) = b1_site_im[2];
		
      poke_real(b22_j) = b2_site_re[2];
      poke_imag(b22_j) = b2_site_im[2];

    } 
    if_dobs3.end();


    poke_real(f0_j) = f_site_re[0];
    poke_imag(f0_j) = f_site_im[0];

    poke_real(f1_j) = f_site_re[1];
    poke_imag(f1_j) = f_site_im[1];

    poke_real(f2_j) = f_site_re[2];
    poke_imag(f2_j) = f_site_im[2];


    auto lv_c0_neg = real_t( c0 < real_t(0.0) ).elem().get_val();
      
    JitIf if_c0_neg2( lv_c0_neg );   // if( c0_negativeP ) 
    {
      poke_imag(f0_j) *= real_t(-1.0);
      poke_real(f1_j) *= real_t(-1.0);
      poke_imag(f2_j) *= real_t(-1.0);
    }
    if_c0_neg2.end();

  }
  c1_lt_0p004.end(); // if (c1 < 4.0e-3 )


  jit_get_function(function);
}



#endif
#endif
