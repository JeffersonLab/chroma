#include "qdp.h"

using namespace QDP;


#ifdef QDP_IS_QDPJIT
#ifndef QDPJIT_IS_QDPJITPTX
#ifndef QDPJIT_IS_QDPJITNVVM

#warning "Using QDP-JIT/LLVM stouting routines"

void function_get_fs_bs_exec(const JitFunction& function, 
			     const LatticeColorMatrix& Q,
			     const LatticeColorMatrix& QQ,
			     multi1d<LatticeComplex>& f,
			     multi1d<LatticeComplex>& b1,
			     multi1d<LatticeComplex>& b2,
			     bool dobs)
{
  AddressLeaf addr_leaf(all);

  addr_leaf.setLit( dobs );

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

  //QDPIO::cerr << "calling getFsBs\n";

  jit_dispatch(function.func().at(0),Layout::sitesOnNode(),getDataLayoutInnerSize(),true,0,addr_leaf);
}



WordREG<REAL> jit_constant( double f )
{
  return WordREG<REAL>(f);
}


void function_get_fs_bs_build( JitFunction& func,
			       const LatticeColorMatrix& Q,
			       const LatticeColorMatrix& QQ,
			       multi1d<LatticeComplex>& f,
			       multi1d<LatticeComplex>& b1,
			       multi1d<LatticeComplex>& b2)
{
  //std::cout << __PRETTY_FUNCTION__ << ": entering\n";

  JitMainLoop loop;

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

  llvm::Value*  r_dobs   = llvm_derefParam( p_dobs );

  IndexDomainVector idx = loop.getIdx();
      
  typename LCMJIT::Subtype_t& Q_j  = Q_jit.elem(JitDeviceLayout::LayoutCoalesced,idx);
  typename LCMJIT::Subtype_t& QQ_j = QQ_jit.elem(JitDeviceLayout::LayoutCoalesced,idx);

  typename LCJIT::Subtype_t& f0_j = f0_jit.elem(JitDeviceLayout::LayoutCoalesced,idx);
  typename LCJIT::Subtype_t& f1_j = f1_jit.elem(JitDeviceLayout::LayoutCoalesced,idx);
  typename LCJIT::Subtype_t& f2_j = f2_jit.elem(JitDeviceLayout::LayoutCoalesced,idx);
  
  typename LCJIT::Subtype_t& b10_j = b10_jit.elem(JitDeviceLayout::LayoutCoalesced,idx);
  typename LCJIT::Subtype_t& b11_j = b11_jit.elem(JitDeviceLayout::LayoutCoalesced,idx);
  typename LCJIT::Subtype_t& b12_j = b12_jit.elem(JitDeviceLayout::LayoutCoalesced,idx);
  typename LCJIT::Subtype_t& b20_j = b20_jit.elem(JitDeviceLayout::LayoutCoalesced,idx);
  typename LCJIT::Subtype_t& b21_j = b21_jit.elem(JitDeviceLayout::LayoutCoalesced,idx);
  typename LCJIT::Subtype_t& b22_j = b22_jit.elem(JitDeviceLayout::LayoutCoalesced,idx);

  llvm::BasicBlock * block_exit = llvm_new_basic_block();
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

    llvm::BasicBlock * block_c1_lt = llvm_new_basic_block();
    llvm::BasicBlock * block_c1_nlt = llvm_new_basic_block();
    llvm_cond_branch( llvm_lt( c1.get_val() , llvm_create_value( 4.0e-3 ) ) , block_c1_lt , block_c1_nlt );

    llvm_set_insert_point( block_c1_lt );
    { //    if( c1 < 4.0e-3  ) 
      f0_j.elem().elem().real() = jit_constant(1.0) - c0 * c0 / jit_constant(720.0);
      f0_j.elem().elem().imag() =  -( c0 / jit_constant(6.0) )*( jit_constant(1.0) -(c1/jit_constant(20.0))*(jit_constant(1.0)-(c1/jit_constant(42.0)))) ;

      f1_j.elem().elem().real() =  c0/jit_constant(24.0)*(jit_constant(1.0)-c1/jit_constant(15.0)*(jit_constant(1.0)-jit_constant(3.0)*c1/jit_constant(112.0))) ;
      f1_j.elem().elem().imag() =  jit_constant(1.0)-c1/jit_constant(6.0)*(jit_constant(1.0)-c1/jit_constant(20.0)*(jit_constant(1.0)-c1/jit_constant(42.0)))-c0*c0/jit_constant(5040.0);
	  
      f2_j.elem().elem().real() = jit_constant(0.5)*(jit_constant(-1.0)+c1/jit_constant(12.0)*(jit_constant(1.0)-c1/jit_constant(30.0)*(jit_constant(1.0)-c1/jit_constant(56.0)))+c0*c0/jit_constant(20160.0));
      f2_j.elem().elem().imag() = jit_constant(0.5)*(c0/jit_constant(60.0)*(jit_constant(1.0)-c1/jit_constant(21.0)*(jit_constant(1.0)-c1/jit_constant(48.0))));


      llvm::BasicBlock * block_dobs_0 = llvm_new_basic_block();
      llvm_cond_branch( r_dobs , block_dobs_0 , block_exit );

      { // dobs
	llvm_set_insert_point( block_dobs_0 );

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

	llvm_branch( block_exit );

      } // Dobs==true

    } // if (c1 < 4.0e-3 )

    llvm_set_insert_point( block_c1_nlt );


    llvm::Value*  c0_negativeP = llvm_lt( c0.get_val() , llvm_create_value(0.0) );
    WordREG<REAL> c0abs = fabs(c0);
    WordREG<REAL> c0max = jit_constant(2.0) * pow( c1 / jit_constant(3.0) , jit_constant(1.5) );
    WordREG<REAL> theta;
    WordREG<REAL> eps = (c0max - c0abs)/c0max;


    // llvm::BasicBlock * block_dobs_0 = llvm_new_basic_block();
    // llvm_cond_branch( r_dobs , block_dobs_0 , block_exit );

    llvm::Value* theta_phi0;
    llvm::Value* theta_phi1;
    llvm::Value* theta_phi2;

    llvm::BasicBlock * block_eps_0 = llvm_new_basic_block();
    llvm::BasicBlock * block_eps_1 = llvm_new_basic_block();
    llvm::BasicBlock * block_eps_2 = llvm_new_basic_block();
    llvm::BasicBlock * block_eps_3 = llvm_new_basic_block();
    llvm::BasicBlock * block_eps_exit = llvm_new_basic_block();

    llvm_cond_branch( llvm_lt( eps.get_val() , llvm_create_value(0.0) ) , block_eps_0 , block_eps_1 );

    llvm_set_insert_point( block_eps_0 ); // epsilon < 0
    theta_phi0 = llvm_create_value(0.0);
    llvm_branch(block_eps_exit);

    llvm_set_insert_point( block_eps_1 ); 
    llvm_cond_branch( llvm_lt( eps.get_val() , llvm_create_value(1.0e-3) ) , block_eps_2 , block_eps_3 );

    llvm_set_insert_point( block_eps_2 ); // epsilon < 1.e-3
    WordREG<REAL> sqtwo = sqrt( jit_constant(2.0) );
    WordREG<REAL> theta_tmp;
    theta_tmp = 
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
    theta_phi1 = theta_tmp.get_val();
    llvm_branch(block_eps_exit);

    llvm_set_insert_point( block_eps_3 ); // else
    WordREG<REAL> theta_tmp2;
    theta_tmp2 = acos( c0abs/c0max );
    theta_phi2 = theta_tmp2.get_val();
    llvm_branch(block_eps_exit);

    llvm_set_insert_point( block_eps_exit ); // exit

    llvm::PHINode* theta_phi = llvm_phi( llvm_type<REAL>::value , 3 );
    theta_phi->addIncoming( theta_phi0 , block_eps_0 );
    theta_phi->addIncoming( theta_phi1 , block_eps_2 );
    theta_phi->addIncoming( theta_phi2 , block_eps_3 );

    theta.setup( theta_phi );


	
    multi1d<WordREG<REAL> > f_site_re(3);
    multi1d<WordREG<REAL> > f_site_im(3);
	  
    multi1d<WordREG<REAL> > b1_site_re_phi0(3);
    multi1d<WordREG<REAL> > b1_site_im_phi0(3);
	      
    multi1d<WordREG<REAL> > b2_site_re_phi0(3);
    multi1d<WordREG<REAL> > b2_site_im_phi0(3);
	  
	  
	  
    WordREG<REAL> u = sqrt(c1/jit_constant(3.0))*cos(theta/jit_constant(3.0));
    WordREG<REAL> w = sqrt(c1)*sin(theta/jit_constant(3.0));
	  
    WordREG<REAL> u_sq = u*u;
    WordREG<REAL> w_sq = w*w;
	  
    WordREG<REAL> xi0,xi1;

    llvm::Value *w_smallP = llvm_lt( (fabs( w )).get_val() , llvm_create_value( 0.05 ) );
    llvm::BasicBlock * block_xi0_exit = llvm_new_basic_block();
    { // xi0
      llvm::BasicBlock * block_xi0_small = llvm_new_basic_block();
      llvm::BasicBlock * block_xi0_not_small = llvm_new_basic_block();
      llvm::Value* xi0_phi0;
      llvm::Value* xi0_phi1;

      llvm_cond_branch( w_smallP , block_xi0_small , block_xi0_not_small );

      llvm_set_insert_point( block_xi0_small ); 
      WordREG<REAL> xi0_tmp0 =
	jit_constant(1.0) - 
	(jit_constant(1.0/6.0)*w_sq*( jit_constant(1.0) - 
				      (jit_constant(1.0/20.)*w_sq*( jit_constant(1.0) - 
								    (jit_constant(1.0/42.0)*w_sq ) ))));
      xi0_phi0 = xi0_tmp0.get_val();

      llvm_branch( block_xi0_exit );

      llvm_set_insert_point( block_xi0_not_small ); 
      WordREG<REAL> xi0_tmp1 = sin(w)/w;
      xi0_phi1 = xi0_tmp1.get_val();
      llvm_branch( block_xi0_exit );

      llvm_set_insert_point( block_xi0_exit ); 

      llvm::PHINode* xi0_phi = llvm_phi( llvm_type<REAL>::value , 2 );
      xi0_phi->addIncoming( xi0_phi0 , block_xi0_small );
      xi0_phi->addIncoming( xi0_phi1 , block_xi0_not_small );
      xi0.setup( xi0_phi );
    } // xi0

    llvm::BasicBlock * block_dobs1 = llvm_new_basic_block();
    llvm::BasicBlock * block_dobs1_exit = llvm_new_basic_block();
    llvm_cond_branch( r_dobs , block_dobs1 , block_dobs1_exit );
    llvm_set_insert_point( block_dobs1 );
    //llvm::BasicBlock * block_xi1_exit = llvm_new_basic_block();
    llvm::BasicBlock * block_xi1_small = llvm_new_basic_block();
    llvm::BasicBlock * block_xi1_not_small = llvm_new_basic_block();
    llvm::Value* xi1_phi0;
    llvm::Value* xi1_phi1;
    { // xi1

      llvm_cond_branch( w_smallP , block_xi1_small , block_xi1_not_small );

      llvm_set_insert_point( block_xi1_small ); 
      WordREG<REAL> xi1_tmp0 = 	    
	jit_constant(-1.0)*
	( jit_constant((REAL)1/(REAL)3) - 
	  jit_constant((REAL)1/(REAL)30)*w_sq*( jit_constant((REAL)1) - 
						jit_constant((REAL)1/(REAL)28)*w_sq*( jit_constant((REAL)1) - 
										      jit_constant((REAL)1/(REAL)54)*w_sq ) ) );
      xi1_phi0 = xi1_tmp0.get_val();
      //llvm_branch( block_xi1_exit );
      llvm_branch( block_dobs1_exit );

      llvm_set_insert_point( block_xi1_not_small ); 
      WordREG<REAL> xi1_tmp1 = cos(w)/w_sq - sin(w)/(w_sq*w);
      xi1_phi1 = xi1_tmp1.get_val();
      //llvm_branch( block_xi1_exit );
      llvm_branch( block_dobs1_exit );

      //llvm_set_insert_point( block_xi1_exit ); 

      //llvm_branch( block_dobs1_exit );
    } // xi1
    llvm_set_insert_point( block_dobs1_exit );

    llvm::PHINode* xi1_phi = llvm_phi( llvm_type<REAL>::value , 3 );
    xi1_phi->addIncoming( xi1_phi0 , block_xi1_small );
    xi1_phi->addIncoming( xi1_phi1 , block_xi1_not_small );
    xi1_phi->addIncoming( llvm_cast( llvm_type<REAL>::value , llvm_create_value(0.0) ) , block_xi0_exit );
    xi1.setup( xi1_phi );

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

    llvm::BasicBlock * block_dobs2 = llvm_new_basic_block();
    llvm::BasicBlock * block_dobs2_exit = llvm_new_basic_block();
    llvm_cond_branch( r_dobs , block_dobs2 , block_dobs2_exit );
    llvm_set_insert_point( block_dobs2 );
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
		  
	  b1_site_re_phi0[j]=( subexp1*r_1_re[j] + subexp2*r_2_re[j] - subexp3*f_site_re[j] )/b_denum;
	  b1_site_im_phi0[j]=( subexp1*r_1_im[j] + subexp2*r_2_im[j] - subexp3*f_site_im[j] )/b_denum;
	}
		
	{ 
	  WordREG<REAL> subexp1 = jit_constant(3.0)*u;
	  WordREG<REAL> subexp2 = jit_constant(24.0)*u;

	  b2_site_re_phi0[j]=( r_1_re[j]- subexp1*r_2_re[j] - subexp2 * f_site_re[j] )/b_denum;
	  b2_site_im_phi0[j]=( r_1_im[j] -subexp1*r_2_im[j] - subexp2 * f_site_im[j] )/b_denum;
	}
      }


      // Now flip the coefficients of the b-s

      llvm::BasicBlock * block_c0neg = llvm_new_basic_block();
      llvm::BasicBlock * block_c0neg_exit = llvm_new_basic_block();
      llvm_cond_branch( c0_negativeP , block_c0neg , block_c0neg_exit );
      llvm_set_insert_point( block_c0neg );

      multi1d<WordREG<REAL> > b1_site_re_phi1(3);
      multi1d<WordREG<REAL> > b1_site_im_phi1(3);
	      
      multi1d<WordREG<REAL> > b2_site_re_phi1(3);
      multi1d<WordREG<REAL> > b2_site_im_phi1(3);

      //if( c0_negativeP ) 
      {
	//b1_site[0] = conj(b1_site[0]);
	b1_site_im_phi1[0] = b1_site_im_phi0[0] * jit_constant(-1.0);
		
	//b1_site[1] = -conj(b1_site[1]);
	b1_site_re_phi1[1] = b1_site_re_phi0[1] * jit_constant(-1.0);
		
	//b1_site[2] = conj(b1_site[2]);
	b1_site_im_phi1[2] = b1_site_im_phi0[2] * jit_constant(-1.0);
		
	//b2_site[0] = -conj(b2_site[0]);
	b2_site_re_phi1[0] = b2_site_re_phi0[0] * jit_constant(-1.0);
		
	//b2_site[1] = conj(b2_site[1]);
	b2_site_im_phi1[1] = b2_site_im_phi0[1] * jit_constant(-1.0);
		
	//b2_site[2] = -conj(b2_site[2]);
	b2_site_re_phi1[2] = b2_site_re_phi0[2] * jit_constant(-1.0);
      }
      llvm_branch( block_c0neg_exit );

      llvm_set_insert_point( block_c0neg_exit );

      //
      // Now, PHI' them together
      //
      multi1d<WordREG<REAL> > f_site_re(3);
      multi1d<WordREG<REAL> > f_site_im(3);
	  
      multi1d<WordREG<REAL> > b1_site_re(3);
      multi1d<WordREG<REAL> > b1_site_im(3);
      
      multi1d<WordREG<REAL> > b2_site_re(3);
      multi1d<WordREG<REAL> > b2_site_im(3);

      b1_site_im[1] = b1_site_im_phi0[1];
      b1_site_re[0] = b1_site_re_phi0[0];
      b1_site_re[2] = b1_site_re_phi0[2];

      b2_site_re[1] = b2_site_re_phi0[1];
      b2_site_im[0] = b2_site_im_phi0[0];
      b2_site_im[2] = b2_site_im_phi0[2];


      qdpPHI( b1_site_im[0] , b1_site_im_phi1[0] , block_c0neg ,
	      b1_site_im_phi0[0] , block_dobs2 );
      qdpPHI( b1_site_re[1] , b1_site_re_phi1[1] , block_c0neg ,
	      b1_site_re_phi0[1] , block_dobs2 );
      qdpPHI( b1_site_im[2] , b1_site_im_phi1[2] , block_c0neg ,
	      b1_site_im_phi0[2] , block_dobs2 );
      qdpPHI( b2_site_re[0] , b2_site_re_phi1[0] , block_c0neg ,
	      b2_site_re_phi0[0] , block_dobs2 );
      qdpPHI( b2_site_im[1] , b2_site_im_phi1[1] , block_c0neg ,
	      b2_site_im_phi0[1] , block_dobs2 );
      qdpPHI( b2_site_re[2] , b2_site_re_phi1[2] , block_c0neg ,
	      b2_site_re_phi0[2] , block_dobs2 );

	      
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

      llvm_branch( block_dobs2_exit );
      //}
    } // end of if (dobs==true)
    //llvm_label( cont_6 );

    llvm_set_insert_point( block_dobs2_exit );


    // Now when everything is done flip signs of the b-s (can't do this before
    // as the unflipped f-s are needed to find the b-s
	  
    // Load back into the lattice sized object
    //      for(int j=0; j < 3; j++) { 
    f0_j.elem().elem().real() = f_site_re[0];
    f0_j.elem().elem().imag() = f_site_im[0];

    f1_j.elem().elem().real() = f_site_re[1];
    f1_j.elem().elem().imag() = f_site_im[1];

    f2_j.elem().elem().real() = f_site_re[2];
    f2_j.elem().elem().imag() = f_site_im[2];
    //}


    llvm::BasicBlock * block_c0neg1 = llvm_new_basic_block();
    llvm::BasicBlock * block_c0neg1_exit = llvm_new_basic_block();
    llvm_cond_branch( c0_negativeP , block_c0neg1 , block_c0neg1_exit );
    llvm_set_insert_point( block_c0neg1 );

    //      if( c0_negativeP ) {
    // f_site[0] = conj(f_site[0]);
    //f_site_im[0] *= jit_constant(-1.0);
    f0_j.elem().elem().imag() *= jit_constant(-1.0);
	    
    //f_site[1] = -conj(f_site[1]);
    //f_site_re[1] *= jit_constant(-1.0);
    f1_j.elem().elem().real() *= jit_constant(-1.0);
	    
    //f_site[2] = conj(f_site[2]);
    //f_site_im[2] *= jit_constant(-1.0);
    f2_j.elem().elem().imag() *= jit_constant(-1.0);

    llvm_branch( block_c0neg1_exit );

    llvm_set_insert_point( block_c0neg1_exit );
	    
    //}

  } // End of if( corner_caseP ) else {}

  llvm_branch( block_exit );
  llvm_set_insert_point( block_exit );


  loop.done();

  func.func().push_back( jit_function_epilogue_get("jit_get_fs_bs.ptx") );
}

#endif
#endif
#endif
