#include "qdp.h"

#ifdef QDP_IS_QDPJIT

using namespace QDP;

void function_hyp_smear_link_exec(JitFunction& function, 
                                  const LatticeColorMatrix& u, 
                                  LatticeColorMatrix& u_hyp,
                                  const Real alpha1,
                                  const Real alpha2,
                                  const Real alpha3,
                                  const int BlkMax,
                                  const Real BlkAccu)
{
  AddressLeaf addr_leaf(all);

  int junk_0 = forEach(u, addr_leaf, NullCombine());
  int junk_1 = forEach(u_hyp, addr_leaf, NullCombine());

  std::vector<QDPCache::ArgKey> ids;

  int lo = 0;
  int hi = Layout::sitesOnNode();
  JitParam jit_lo( QDP_get_global_cache().addJitParamInt( lo ) );
  JitParam jit_hi( QDP_get_global_cache().addJitParamInt( hi ) );
  ids.push_back( jit_lo.get_id() );
  ids.push_back( jit_hi.get_id() );

  JitParam jit_alpha1( QDP_get_global_cache().addJitParamDouble( toDouble(alpha1) ) );
  JitParam jit_alpha2( QDP_get_global_cache().addJitParamDouble( toDouble(alpha2) ) );
  JitParam jit_alpha3( QDP_get_global_cache().addJitParamDouble( toDouble(alpha3) ) );
  JitParam jit_BlkAccu( QDP_get_global_cache().addJitParamDouble( toDouble(BlkAccu) ) );
  JitParam jit_BlkMax( QDP_get_global_cache().addJitParamInt( BlkMax ) );
  
  ids.push_back( jit_alpha1.get_id() );
  ids.push_back( jit_alpha2.get_id() );
  ids.push_back( jit_alpha3.get_id() );
  ids.push_back( jit_BlkAccu.get_id() );
  ids.push_back( jit_BlkMax.get_id() );

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

void function_hyp_smear_link_build(JitFunction& function,
                                   const LatticeColorMatrix& u, 
                                   LatticeColorMatrix& u_hyp,
                                   const Real alpha1,
                                   const Real alpha2,
                                   const Real alpha3,
                                   const int BlkMax,
                                   const Real BlkAccu)
{
  if (ptx_db::db_enabled)
    {
      llvm_ptx_db( function , __PRETTY_FUNCTION__ );
      if (!function.empty())
	return;
    }

  llvm_start_new_function("hyp_smear_link",__PRETTY_FUNCTION__);

  ParamRef  p_lo     = llvm_add_param<int>();
  ParamRef  p_hi     = llvm_add_param<int>();
  ParamRef  p_alpha1 = llvm_add_param<REAL>();
  ParamRef  p_alpha2 = llvm_add_param<REAL>();
  ParamRef  p_alpha3 = llvm_add_param<REAL>();
  ParamRef  p_BlkAccu = llvm_add_param<REAL>();
  ParamRef  p_BlkMax = llvm_add_param<int>();

  ParamLeaf param_leaf;
  
  typedef typename LeafFunctor<LatticeColorMatrix, ParamLeaf>::Type_t  LCMJIT;
  typedef typename LeafFunctor<LatticeComplex    , ParamLeaf>::Type_t  LCJIT;

  LCMJIT u_jit(forEach(u, param_leaf, TreeCombine()));
  LCMJIT u_hyp_jit(forEach(u_hyp, param_leaf, TreeCombine()));

  llvm::Value*  r_lo     = llvm_derefParam( p_lo );
  llvm::Value*  r_hi     = llvm_derefParam( p_hi );
  llvm::Value*  r_alpha1 = llvm_derefParam( p_alpha1 );
  llvm::Value*  r_alpha2 = llvm_derefParam( p_alpha2 );
  llvm::Value*  r_alpha3 = llvm_derefParam( p_alpha3 );
  llvm::Value*  r_BlkAccu= llvm_derefParam( p_BlkAccu );
  llvm::Value*  r_BlkMax = llvm_derefParam( p_BlkMax );
      
  llvm::Value*  r_idx = llvm_thread_idx();
      
  llvm_cond_exit( llvm_ge( r_idx , r_hi ) );


  auto u_j  = u_jit.elem(JitDeviceLayout::Coalesced,r_idx);
  auto u_hyp_j = u_hyp_jit.elem(JitDeviceLayout::Coalesced,r_idx);

  jit_get_function(function);
}
#endif
