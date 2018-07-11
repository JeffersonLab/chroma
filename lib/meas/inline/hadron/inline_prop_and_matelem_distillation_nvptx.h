#ifndef inline_prop_and_matelem_distillation_nvptx_h
#define inline_prop_and_matelem_distillation_nvptx_h

namespace QDP
{


  // T1 input
  // T2 output
  template< class T1 , class T2 , JitDeviceLayout input_layout >
  CUfunction 
  function_multi_sum_convert_build()
  {
    if (ptx_db::db_enabled) {
      CUfunction func = llvm_ptx_db( __PRETTY_FUNCTION__ );
      if (func)
	return func;
    }

    llvm_start_new_function();

    ParamRef p_lo        = llvm_add_param<int>();
    ParamRef p_hi        = llvm_add_param<int>();
    ParamRef p_inputnum  = llvm_add_param<int>();

    typedef typename WordType<T1>::Type_t T1WT;
    typedef typename WordType<T2>::Type_t T2WT;

    ParamRef p_idata      = llvm_add_param< T1WT** >();  // Input  array
    ParamRef p_odata      = llvm_add_param< T2WT** >();  // output array

    OLatticeJIT<typename JITType<T1>::Type_t> idata(  p_idata );   // want coal   access later
    OLatticeJIT<typename JITType<T2>::Type_t> odata(  p_odata );   // want scalar access later

    llvm_derefParam( p_lo ); // r_lo
    llvm::Value* r_hi           = llvm_derefParam( p_hi );
    llvm::Value* r_inputnum     = llvm_derefParam( p_inputnum );

    llvm::Value* r_shared = llvm_get_shared_ptr( llvm_type<T2WT>::value );

    typedef typename JITType<T2>::Type_t T2JIT;

    llvm::Value* r_idx = llvm_thread_idx();

    llvm::Value* r_block_idx  = llvm_call_special_ctaidx();
    llvm::Value* r_tidx       = llvm_call_special_tidx();
    llvm::Value* r_ntidx       = llvm_call_special_ntidx(); // needed later

    IndexDomainVector args;
    args.push_back( make_pair( Layout::sitesOnNode() , r_tidx ) );  // sitesOnNode irrelevant since Scalar access later
    T2JIT sdata_jit;
    sdata_jit.setup( r_shared , JitDeviceLayout::Scalar , args );
    zero_rep( sdata_jit );

    llvm_cond_exit( llvm_ge( r_idx , r_hi ) );

    llvm::BasicBlock * entry_block = llvm_get_insert_block();

    //
    // big loop over array of input vectors,   0 <= r_input < r_inputnum
    //
    llvm::BasicBlock * block_input_loop_start = llvm_new_basic_block();
    llvm::BasicBlock * block_input_loop_inc = llvm_new_basic_block();
    llvm::BasicBlock * block_input_loop_exit = llvm_new_basic_block();
    llvm::Value* r_input_inc;

    llvm_branch( block_input_loop_start );

    llvm_set_insert_point( block_input_loop_start );

    llvm::PHINode * r_input = llvm_phi( llvm_type<int>::value , 2 );
    r_input->addIncoming( llvm_create_value(0) , entry_block );

    llvm_cond_branch( llvm_ge( r_input , r_inputnum ) , block_input_loop_exit , block_input_loop_inc );
    {
      llvm_set_insert_point(block_input_loop_inc);
    
      typename REGType< typename JITType<T1>::Type_t >::Type_t reg_idata_elem;   // this is stupid
      reg_idata_elem.setup( idata.elem( input_layout , r_idx , r_input ) );

      sdata_jit = reg_idata_elem; // This should do the precision conversion (SP->DP)

      llvm_bar_sync();

      //llvm::Value* val_ntid = llvm_call_special_ntidx();

      //
      // Find next power of 2 loop
      //
      llvm::BasicBlock * block_power_loop_start = llvm_new_basic_block();
      llvm::BasicBlock * block_power_loop_inc = llvm_new_basic_block();
      llvm::BasicBlock * block_power_loop_exit = llvm_new_basic_block();
      llvm::Value* r_pow_phi;

      llvm_branch( block_power_loop_start );

      llvm_set_insert_point( block_power_loop_start );

      llvm::PHINode * r_pow = llvm_phi( llvm_type<int>::value , 2 );
      r_pow->addIncoming( llvm_create_value(1) , block_input_loop_inc );   // entry_block, block_input_loop_exit

      llvm_cond_branch( llvm_ge( r_pow , r_ntidx ) , block_power_loop_exit , block_power_loop_inc );
      {
	llvm_set_insert_point(block_power_loop_inc);
	r_pow_phi = llvm_shl( r_pow , llvm_create_value(1) );
	r_pow->addIncoming( r_pow_phi , block_power_loop_inc );
	llvm_branch( block_power_loop_start );
      }

      llvm_set_insert_point(block_power_loop_exit);

      llvm::Value* r_pow_shr1 = llvm_shr( r_pow , llvm_create_value(1) );

      //
      // Shared memory reduction loop
      //
      llvm::BasicBlock * block_red_loop_start = llvm_new_basic_block();
      llvm::BasicBlock * block_red_loop_start_1 = llvm_new_basic_block();
      llvm::BasicBlock * block_red_loop_start_2 = llvm_new_basic_block();
      llvm::BasicBlock * block_red_loop_add = llvm_new_basic_block();
      llvm::BasicBlock * block_red_loop_sync = llvm_new_basic_block();
      llvm::BasicBlock * block_red_loop_end = llvm_new_basic_block();

      llvm_branch( block_red_loop_start );
      llvm_set_insert_point(block_red_loop_start);
    
      llvm::PHINode * r_red_pow = llvm_phi( llvm_type<int>::value , 2 );    
      r_red_pow->addIncoming( r_pow_shr1 , block_power_loop_exit );
      llvm_cond_branch( llvm_le( r_red_pow , llvm_create_value(0) ) , block_red_loop_end , block_red_loop_start_1 );

      llvm_set_insert_point(block_red_loop_start_1);

      llvm_cond_branch( llvm_ge( r_tidx , r_red_pow ) , block_red_loop_sync , block_red_loop_start_2 );

      llvm_set_insert_point(block_red_loop_start_2);

      llvm::Value * v = llvm_add( r_red_pow , r_tidx );
      llvm_cond_branch( llvm_ge( v , r_ntidx ) , block_red_loop_sync , block_red_loop_add );

      llvm_set_insert_point(block_red_loop_add);


      IndexDomainVector args_new;
      args_new.push_back( make_pair( Layout::sitesOnNode() , 
				     llvm_add( r_tidx , r_red_pow ) ) );  // sitesOnNode irrelevant since Scalar access later

      typename JITType<T2>::Type_t sdata_jit_plus;
      sdata_jit_plus.setup( r_shared , JitDeviceLayout::Scalar , args_new );

      typename REGType< typename JITType<T2>::Type_t >::Type_t sdata_reg_plus;    // 
      sdata_reg_plus.setup( sdata_jit_plus );

      sdata_jit += sdata_reg_plus;


      llvm_branch( block_red_loop_sync );

      llvm_set_insert_point(block_red_loop_sync);
      llvm_bar_sync();
      llvm::Value* pow_1 = llvm_shr( r_red_pow , llvm_create_value(1) );
      r_red_pow->addIncoming( pow_1 , block_red_loop_sync );

      llvm_branch( block_red_loop_start );

      llvm_set_insert_point(block_red_loop_end);


      llvm::BasicBlock * block_store_global = llvm_new_basic_block();
      llvm::BasicBlock * block_not_store_global = llvm_new_basic_block();
      llvm_cond_branch( llvm_eq( r_tidx , llvm_create_value(0) ) , 
			block_store_global , 
			block_not_store_global );
      llvm_set_insert_point(block_store_global);
      typename REGType< typename JITType<T2>::Type_t >::Type_t sdata_reg;   // this is stupid
      sdata_reg.setup( sdata_jit );
      odata.elem( JitDeviceLayout::Scalar , r_block_idx , r_input ) = sdata_reg;

      llvm_branch( block_not_store_global );
      llvm_set_insert_point(block_not_store_global);
      
      r_input_inc = llvm_add( r_input , llvm_create_value(1) );
      r_input->addIncoming( r_input_inc , block_not_store_global );  // block_input_loop_inc
    
      llvm_branch( block_input_loop_start );
    }

    llvm_set_insert_point(block_input_loop_exit);
    
    //llvm_module_dump();
    
    return jit_function_epilogue_get_cuf("multi_sum.ptx" , __PRETTY_FUNCTION__ );
  }


  template < class T1 , class T2 , JitDeviceLayout input_layout >
  void qdp_jit_multi_reduce_convert(int size, 
				    int threads, 
				    int blocks, 
				    int shared_mem_usage,
				    multi1d<T1*>& d_idata, 
				    multi1d<T2*>& d_odata )
  {
    static CUfunction function;

    // Build the function
    if (function == NULL)
      {
	//std::cout << __PRETTY_FUNCTION__ << ": does not exist - will build\n";
	function = function_multi_sum_convert_build<T1,T2,input_layout>();
	//std::cout << __PRETTY_FUNCTION__ << ": did not exist - finished building\n";
      }
    else
      {
	//std::cout << __PRETTY_FUNCTION__ << ": is already built\n";
      }

    // Execute the function
    function_multi_sum_convert_exec(function, size, threads, blocks, shared_mem_usage, 
     				    d_idata, d_odata );
  }

  template < class T1 , class T2 >
  void
  function_multi_sum_convert_exec( CUfunction function, 
				   int size, int threads, int blocks, int shared_mem_usage,
				   multi1d<T1*>& d_idata, multi1d<T2*>& d_odata )
  {
    //std::cout << __PRETTY_FUNCTION__ << ": entering\n";

    assert( d_idata.size() == d_odata.size() );

    const unsigned N = d_idata.size();
    
    int lo = 0;
    int hi = size;
    
    int inputnum = d_idata.size();

    std::vector<void*> addr;

    addr.push_back( &lo );
    addr.push_back( &hi );
    addr.push_back( &inputnum );

    void* ptr_i;
    void* ptr_o;

    if (!QDP_get_global_cache().allocate_device_static( (void**)&(ptr_i) , N*sizeof(T1*) ))
      QDP_error_exit( "multi_sum buffer: buffer no memory, exit");
    if (!QDP_get_global_cache().allocate_device_static( (void**)&(ptr_o) , N*sizeof(T1*) ))
      QDP_error_exit( "multi_sum buffer: buffer no memory, exit");

    CudaMemcpyH2D( ptr_i , d_idata.slice() , N*sizeof(T1*) );
    CudaMemcpyH2D( ptr_o , d_odata.slice() , N*sizeof(T1*) );

    addr.push_back( &ptr_i );
    addr.push_back( &ptr_o );
    
    kernel_geom_t now = getGeom( hi-lo , threads );

    //QDP_info("launing block=(%d,1,1)  grid=(%d,%d,1)",threads,now.Nblock_x,now.Nblock_y);

    CudaLaunchKernel(function,   now.Nblock_x,now.Nblock_y,1,    threads,1,1,    shared_mem_usage, 0, &addr[0] , 0);

    QDP_get_global_cache().free_device_static( ptr_i );
    QDP_get_global_cache().free_device_static( ptr_o );
  }


  template < class T1 , class T2 , class T3 >
  void
  function_multi_localInnerProduct_sum_convert_exec( CUfunction function, 
						     int size, int threads, int blocks, int shared_mem_usage,
						     multi1d<T1*>& d_idata, multi1d<T2*>& d_odata , T3* d_vdata ,
						     const Subset& s )
  {
    //std::cout << __PRETTY_FUNCTION__ << ": entering\n";

    void * site_table = QDP_get_global_cache().getDevicePtr( s.getIdSiteTable() );

    assert( d_idata.size() == d_odata.size() );

    const unsigned N = d_idata.size();
    
    int lo = 0;
    int hi = size;
    
    int inputnum = d_idata.size();

    std::vector<void*> addr;

    addr.push_back( &lo );
    addr.push_back( &hi );
    addr.push_back( &inputnum );

    void* ptr_i;
    void* ptr_o;

    if (!QDP_get_global_cache().allocate_device_static( (void**)&(ptr_i) , N*sizeof(T1*) ))
      QDP_error_exit( "multi_sum buffer: buffer no memory, exit");
    if (!QDP_get_global_cache().allocate_device_static( (void**)&(ptr_o) , N*sizeof(T1*) ))
      QDP_error_exit( "multi_sum buffer: buffer no memory, exit");

    CudaMemcpyH2D( ptr_i , d_idata.slice() , N*sizeof(T1*) );
    CudaMemcpyH2D( ptr_o , d_odata.slice() , N*sizeof(T1*) );

    addr.push_back( &ptr_i );
    addr.push_back( &ptr_o );
    addr.push_back( &d_vdata );
    addr.push_back( &site_table );
    
    kernel_geom_t now = getGeom( hi-lo , threads );

    //QDP_info("launing block=(%d,1,1)  grid=(%d,%d,1)",threads,now.Nblock_x,now.Nblock_y);

    CudaLaunchKernel(function,   now.Nblock_x,now.Nblock_y,1,    threads,1,1,    shared_mem_usage, 0, &addr[0] , 0);

    QDP_get_global_cache().free_device_static( ptr_i );
    QDP_get_global_cache().free_device_static( ptr_o );
  }


  // T1 input
  // T2 output
  // T3 QDPType for localInnerProduct = vector
  template< class T1 , class T2 , class T3 , JitDeviceLayout input_layout >
  CUfunction 
  function_multi_localInnerProduct_sum_convert_build()
  {
    if (ptx_db::db_enabled) {
      CUfunction func = llvm_ptx_db( __PRETTY_FUNCTION__ );
      if (func)
	return func;
    }

    llvm_start_new_function();

    ParamRef p_lo        = llvm_add_param<int>();
    ParamRef p_hi        = llvm_add_param<int>();
    ParamRef p_inputnum  = llvm_add_param<int>();

    typedef typename WordType<T1>::Type_t T1WT;
    typedef typename WordType<T2>::Type_t T2WT;
    typedef typename WordType<T3>::Type_t T3WT;

    ParamRef p_idata      = llvm_add_param< T1WT** >();  // Input  array
    ParamRef p_odata      = llvm_add_param< T2WT** >();  // output array
    ParamRef p_vdata      = llvm_add_param< T3WT* >();  // Vector array

    ParamRef p_site_table   = llvm_add_param<int*>();      // subset sitetable


    OLatticeJIT<typename JITType<T1>::Type_t> idata(  p_idata );   // want scalar access later
    OLatticeJIT<typename JITType<T2>::Type_t> odata(  p_odata );   // want scalar access later
    OLatticeJIT<typename JITType<T3>::Type_t> vdata(  p_vdata );   // want coal   access later

    llvm_derefParam( p_lo ); // r_lo
    llvm::Value* r_hi           = llvm_derefParam( p_hi );
    llvm::Value* r_inputnum     = llvm_derefParam( p_inputnum );

    llvm::Value* r_shared = llvm_get_shared_ptr( llvm_type<T2WT>::value );

    typedef typename JITType<T2>::Type_t T2JIT;

    llvm::Value* r_idx = llvm_thread_idx();
    llvm::Value* r_idx_perm = llvm_array_type_indirection( p_site_table , r_idx );

    llvm::Value* r_block_idx  = llvm_call_special_ctaidx();
    llvm::Value* r_tidx       = llvm_call_special_tidx();
    llvm::Value* r_ntidx       = llvm_call_special_ntidx(); // needed later

    IndexDomainVector args;
    args.push_back( make_pair( Layout::sitesOnNode() , r_tidx ) );  // sitesOnNode irrelevant since Scalar access later
    T2JIT sdata_jit;
    sdata_jit.setup( r_shared , JitDeviceLayout::Scalar , args );
    zero_rep( sdata_jit );

    llvm_cond_exit( llvm_ge( r_idx , r_hi ) );

    llvm::BasicBlock * entry_block = llvm_get_insert_block();

    //
    // big loop over array of input vectors,   0 <= r_input < r_inputnum
    //
    llvm::BasicBlock * block_input_loop_start = llvm_new_basic_block();
    llvm::BasicBlock * block_input_loop_inc = llvm_new_basic_block();
    llvm::BasicBlock * block_input_loop_exit = llvm_new_basic_block();
    llvm::Value* r_input_inc;

    llvm_branch( block_input_loop_start );

    llvm_set_insert_point( block_input_loop_start );

    llvm::PHINode * r_input = llvm_phi( llvm_type<int>::value , 2 );
    r_input->addIncoming( llvm_create_value(0) , entry_block );

    llvm_cond_branch( llvm_ge( r_input , r_inputnum ) , block_input_loop_exit , block_input_loop_inc );
    {
      llvm_set_insert_point(block_input_loop_inc);

      ParamLeaf param_leaf;
      FnLocalInnerProduct op;
      auto op_jit = AddOpParam<FnLocalInnerProduct,ParamLeaf>::apply(op,param_leaf);

      typename REGType< typename JITType<T1>::Type_t >::Type_t reg_idata_elem;
      reg_idata_elem.setup( idata.elem( input_layout , r_idx , r_input ) );

      typename REGType< typename JITType<T3>::Type_t >::Type_t reg_vdata_elem;
      reg_vdata_elem.setup( vdata.elem( JitDeviceLayout::Coalesced , r_idx_perm ) );

      sdata_jit = op_jit( reg_idata_elem , reg_vdata_elem ); // This should do the precision conversion (SP->DP)

      llvm_bar_sync();


      //llvm::Value* val_ntid = llvm_call_special_ntidx();

      //
      // Find next power of 2 loop
      //
      llvm::BasicBlock * block_power_loop_start = llvm_new_basic_block();
      llvm::BasicBlock * block_power_loop_inc = llvm_new_basic_block();
      llvm::BasicBlock * block_power_loop_exit = llvm_new_basic_block();
      llvm::Value* r_pow_phi;

      llvm_branch( block_power_loop_start );

      llvm_set_insert_point( block_power_loop_start );

      llvm::PHINode * r_pow = llvm_phi( llvm_type<int>::value , 2 );
      r_pow->addIncoming( llvm_create_value(1) , block_input_loop_inc );   // entry_block, block_input_loop_exit

      llvm_cond_branch( llvm_ge( r_pow , r_ntidx ) , block_power_loop_exit , block_power_loop_inc );
      {
	llvm_set_insert_point(block_power_loop_inc);
	r_pow_phi = llvm_shl( r_pow , llvm_create_value(1) );
	r_pow->addIncoming( r_pow_phi , block_power_loop_inc );
	llvm_branch( block_power_loop_start );
      }

      llvm_set_insert_point(block_power_loop_exit);

      llvm::Value* r_pow_shr1 = llvm_shr( r_pow , llvm_create_value(1) );

      //
      // Shared memory reduction loop
      //
      llvm::BasicBlock * block_red_loop_start = llvm_new_basic_block();
      llvm::BasicBlock * block_red_loop_start_1 = llvm_new_basic_block();
      llvm::BasicBlock * block_red_loop_start_2 = llvm_new_basic_block();
      llvm::BasicBlock * block_red_loop_add = llvm_new_basic_block();
      llvm::BasicBlock * block_red_loop_sync = llvm_new_basic_block();
      llvm::BasicBlock * block_red_loop_end = llvm_new_basic_block();

      llvm_branch( block_red_loop_start );
      llvm_set_insert_point(block_red_loop_start);
    
      llvm::PHINode * r_red_pow = llvm_phi( llvm_type<int>::value , 2 );    
      r_red_pow->addIncoming( r_pow_shr1 , block_power_loop_exit );
      llvm_cond_branch( llvm_le( r_red_pow , llvm_create_value(0) ) , block_red_loop_end , block_red_loop_start_1 );

      llvm_set_insert_point(block_red_loop_start_1);

      llvm_cond_branch( llvm_ge( r_tidx , r_red_pow ) , block_red_loop_sync , block_red_loop_start_2 );

      llvm_set_insert_point(block_red_loop_start_2);

      llvm::Value * v = llvm_add( r_red_pow , r_tidx );
      llvm_cond_branch( llvm_ge( v , r_ntidx ) , block_red_loop_sync , block_red_loop_add );

      llvm_set_insert_point(block_red_loop_add);


      IndexDomainVector args_new;
      args_new.push_back( make_pair( Layout::sitesOnNode() , 
				     llvm_add( r_tidx , r_red_pow ) ) );  // sitesOnNode irrelevant since Scalar access later

      typename JITType<T2>::Type_t sdata_jit_plus;
      sdata_jit_plus.setup( r_shared , JitDeviceLayout::Scalar , args_new );

      typename REGType< typename JITType<T2>::Type_t >::Type_t sdata_reg_plus;    // 
      sdata_reg_plus.setup( sdata_jit_plus );

      sdata_jit += sdata_reg_plus;


      llvm_branch( block_red_loop_sync );

      llvm_set_insert_point(block_red_loop_sync);
      llvm_bar_sync();
      llvm::Value* pow_1 = llvm_shr( r_red_pow , llvm_create_value(1) );
      r_red_pow->addIncoming( pow_1 , block_red_loop_sync );

      llvm_branch( block_red_loop_start );

      llvm_set_insert_point(block_red_loop_end);


      llvm::BasicBlock * block_store_global = llvm_new_basic_block();
      llvm::BasicBlock * block_not_store_global = llvm_new_basic_block();
      llvm_cond_branch( llvm_eq( r_tidx , llvm_create_value(0) ) , 
			block_store_global , 
			block_not_store_global );
      llvm_set_insert_point(block_store_global);
      typename REGType< typename JITType<T2>::Type_t >::Type_t sdata_reg;   // this is stupid
      sdata_reg.setup( sdata_jit );
      odata.elem( JitDeviceLayout::Scalar , r_block_idx , r_input ) = sdata_reg;

      llvm_branch( block_not_store_global );
      llvm_set_insert_point(block_not_store_global);
      
      r_input_inc = llvm_add( r_input , llvm_create_value(1) );
      r_input->addIncoming( r_input_inc , block_not_store_global );  // block_input_loop_inc
    
      llvm_branch( block_input_loop_start );
    }

    llvm_set_insert_point(block_input_loop_exit);
    
    //llvm_module_dump();
    
    return jit_function_epilogue_get_cuf("multi_sum.ptx" , __PRETTY_FUNCTION__ );
  }





  template < class T1 , class T2 , class T3, JitDeviceLayout input_layout >
  void qdp_jit_multi_localInnerProduct_reduce_convert(int size, 
						      int threads, 
						      int blocks, 
						      int shared_mem_usage,
						      multi1d<T1*>& d_idata, 
						      multi1d<T2*>& d_odata,
						      T3* d_vdata,
						      const Subset& s)
  {
    static CUfunction function;

    // Build the function
    if (function == NULL)
      {
	//std::cout << __PRETTY_FUNCTION__ << ": does not exist - will build\n";
	function = function_multi_localInnerProduct_sum_convert_build<T1,T2,T3,input_layout>();
	//std::cout << __PRETTY_FUNCTION__ << ": did not exist - finished building\n";
      }
    else
      {
	//std::cout << __PRETTY_FUNCTION__ << ": is already built\n";
      }

    // Execute the function
    function_multi_localInnerProduct_sum_convert_exec(function, size, threads, blocks, shared_mem_usage, 
						      d_idata, d_odata, d_vdata , s );
  }

  


  
  

  template<class T, class T3>
  void
  multi_innerProduct( multi1d< ComplexD* >& ret , const multi1d< OSubLattice<T>* >& ms1 , const OLattice<T3>& v1 )
  {
    const int N = ms1.size();

    //QDPIO::cout << "multi_innerProduct (GPU) with N = " << N << "\n";

    assert( N > 0 );
    assert( ret.size() == N );

    for (int i = 0 ; i < N ; ++i )
      {
	if (!ms1[i]->getOwnsMemory())
	  QDP_error_exit("sum with subtype view called");
	assert( ms1[0]->subset().numSiteTable() == ms1[i]->subset().numSiteTable() );
      }

    //typedef typename UnaryReturn<OLattice<T>, FnSum>::Type_t Ret_t;
    //multi1d<Ret_t> ret(N);

    typedef typename BinaryReturn<T,T3,FnLocalInnerProduct>::Type_t       TT3;
    typedef typename UnaryReturn<OLattice<TT3>, FnSum>::Type_t::SubType_t T2;
    
    // these are for the GPU summation
    multi1d<T2*> out_dev(N);
    multi1d<T2*> in_dev(N);

    multi1d<T*>  dev_ptrs(N); // lattice input
    multi1d<T2*> dev_sum(N);  // final   output
    T3* dev_vec;              // lattice input vector for contraction
    
    dev_vec = (T3*)QDP_get_global_cache().getDevicePtr( v1.getId() );
    for (int i = 0 ; i < N ; ++i ) {
      dev_ptrs[i] = (T*)QDP_get_global_cache().getDevicePtr( ms1[i]->getId() );
      dev_sum[i]  = (T2*)QDP_get_global_cache().getDevicePtr( ret[i]->getId() );
      //QDPIO::cout << "input lattice = " << (size_t)dev_ptrs[i] << "\n";
    }

    unsigned actsize=ms1[0]->subset().numSiteTable();
    bool first=true;
    while (1) {

      unsigned numThreads = DeviceParams::Instance().getMaxBlockX();
      while ((numThreads*sizeof(T2) > DeviceParams::Instance().getMaxSMem()) || (numThreads > actsize)) {
	numThreads >>= 1;
      }
      unsigned numBlocks=(int)ceil(float(actsize)/numThreads);
    
      if (numBlocks > DeviceParams::Instance().getMaxGridX()) {
	QDP_error_exit( "sum(Lat,subset) numBlocks(%d) > maxGridX(%d)",numBlocks,(int)DeviceParams::Instance().getMaxGridX());
      }

      int shared_mem_usage = numThreads*sizeof(T2);
      //QDP_info("multi_sum(SubLat): using %d threads per block, %d blocks, shared mem=%d" , numThreads , numBlocks , shared_mem_usage );

      if (first) {
	//QDPIO::cout << "allocating device memory N times\n";
	for (int i = 0 ; i < N ; ++i ) {      
	  if (!QDP_get_global_cache().allocate_device_static( (void**)&(out_dev[i]) , numBlocks*sizeof(T2) ))
	    QDP_error_exit( "sum(lat,subset) reduction buffer: 1st buffer no memory, exit");
	  if (!QDP_get_global_cache().allocate_device_static( (void**)&(in_dev[i]) , numBlocks*sizeof(T2) ))
	    QDP_error_exit( "sum(lat,subset) reduction buffer: 2nd buffer no memory, exit");
	  //QDPIO::cout << "out = " << (size_t)out_dev[i] << "    in = " << (size_t)in_dev[i] << "\n";
	}
      }


      if (numBlocks == 1)
	{
	  if (first)
	    {
	      qdp_jit_multi_localInnerProduct_reduce_convert<T,T2,T3,JitDeviceLayout::Scalar>(actsize, numThreads, numBlocks, shared_mem_usage ,  // ok: Scalar
											      dev_ptrs,
											      dev_sum ,
											      dev_vec ,
											      ms1[0]->subset() );
	    }
	  else
	    {
	      qdp_jit_multi_reduce_convert<T2,T2,JitDeviceLayout::Scalar>( actsize , numThreads , numBlocks, shared_mem_usage , 
									   in_dev , dev_sum );
	    }
	}
      else
	{
	  if (first)
	    {
	      qdp_jit_multi_localInnerProduct_reduce_convert<T,T2,T3,JitDeviceLayout::Scalar>(actsize, numThreads, numBlocks, shared_mem_usage,       // ok: Scalar
											      dev_ptrs,
											      out_dev ,
											      dev_vec ,
											      ms1[0]->subset() );
	    }
	  else
	    {
	      qdp_jit_multi_reduce_convert<T2,T2,JitDeviceLayout::Scalar>( actsize , numThreads , numBlocks, shared_mem_usage , 
									   in_dev , out_dev );
	    }
	}

    
  
      first =false;
    
      if (numBlocks==1) 
	break;

      actsize=numBlocks;
    
      multi1d<T2*> tmp = in_dev;
      for (int i = 0 ; i < N ; ++i ) {      
	in_dev[i] = out_dev[i];
	out_dev[i] = tmp[i];
      }
    }

    //QDPIO::cout << "freeing device memory N times\n";
    for (int i = 0 ; i < N ; ++i ) {      
      QDP_get_global_cache().free_device_static( in_dev[i] );
      QDP_get_global_cache().free_device_static( out_dev[i] );
    }
  
    //QDPIO::cout << "global sum N times\n";
    for (int i = 0 ; i < N ; ++i ) {
      QDPInternal::globalSum(*ret[i]);
      //QDPIO::cout << *ret[i] << "\n";
    }
    
    //return d;
  }


  
} // QDP


#endif
