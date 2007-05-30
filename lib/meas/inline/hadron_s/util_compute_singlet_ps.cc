//
// Wrapper code to compute the connected part of the 
// pseudocalar singlet meson (4 link operator)
//
//

#include "handle.h"
//#include "actions/ferm/fermbcs/fermbcs.h"
#include "actions/ferm/fermacts/fermacts_s.h"
#include "meas/hadron/hadron_s.h"
#include "meas/smear/fuzz_smear.h"
#include "meas/sources/srcfil.h"
#include "util/ferm/transf.h"
#include "meas/hadron/pion_sing_s.h"


#include "util_compute_quark_prop_s.h"

namespace Chroma { 

  enum func_flag{
    QPROP_FUZZ,
    QPROP_NO_FUZZ,
    SINGLETS
  };
  typedef func_flag func_flag_type;


  int check_qprop_source_compatability(stag_src_type type_of_src,
				       bool gauge_shift,
				       bool sym_shift, func_flag_type fflag);

  /***************************************************************************/

  int compute_singlet_ps(LatticeStaggeredFermion & psi,
			 LatticeStaggeredPropagator quark_propagator,
			 stag_src_type type_of_src,
			 bool gauge_shift,
			 bool sym_shift,
			 const multi1d<LatticeColorMatrix> & u ,
			 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
			 XMLWriter & xml_out,
			 Real RsdCG, Real Mass, 
			 int j_decay, int t_source, int t_length){

    LatticeStaggeredFermion q_source ;
    LatticeStaggeredFermion q_source_fuzz ;
    int ncg_had = 0 ;

    LatticeStaggeredPropagator quark_propagator_4link ;

    q_source = zero ;

    // safety checks
    check_qprop_source_compatability(type_of_src, gauge_shift, sym_shift,
				     SINGLETS);

    if( gauge_shift ){

      if((type_of_src == LOCAL_SRC) || (type_of_src == FUZZED_SRC)){
	QDPIO::cerr << "Conflicting source and shift types in " <<
	  "util_compute_quark_prop_s.cc" <<endl;
	exit(0);
      }


    }else{

      if (sym_shift){
	QDPIO::cerr << "Conflicting source and shift types in " <<
	  "util_compute_singlet_s.cc" <<endl;
	QDPIO::cerr << "symmetric shifts not implemented for LLL " <<
	  "non-gauge-invariant sources" << endl;

	if(type_of_src == LOCAL_SRC){
	  QDPIO::cerr << "You asked for LOCAL_SRC, sym_shift" << endl;
	}
	if(type_of_src == FUZZED_SRC){
	  QDPIO::cerr << "You asked for FUZZED_SRC, sym_shift" << endl;
	}

	exit(0);
      }
    
    
      if(type_of_src == GAUGE_INVAR_LOCAL_SOURCE){
	QDPIO::cout << "Conflicting source and shift types in " <<
	  "util_compute_quark_prop_s.cc" <<endl;
	exit(0);
      }
    }


    /*** generate the 4-link quark propagator ****/
    push(xml_out,"Computation_4link_pseudoscalar");


    for(int color_source = 0; color_source < Nc; ++color_source) {

      if( type_of_src == LOCAL_SRC ){

	q_source = zero ;
	multi1d<int> coord(Nd);
	coord[0]=1; coord[1] = 1; coord[2] = 1; coord[3] = 1;
	srcfil(q_source, coord,color_source ) ;
      } else {
	if( type_of_src == GAUGE_INVAR_LOCAL_SOURCE  ) {

	  q_source = zero ;
	  multi1d<int> coord(Nd);
	  
	  // start with local source 
	  coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
	  srcfil(q_source, coord,color_source ) ;


	  
	  // now do the shift
	  coord[0]=1; coord[1] = 1; coord[2] = 1; coord[3] = 1;
	  q_source_fuzz = q_source  ;
	  q_source = shiftDeltaPropCov(coord,q_source_fuzz,u,
				       sym_shift);

	}
      }

      // we fuzz the local quark propagator so it does not
      // make sense to fuzz the shited source quark propagator


      // Compute the propagator for given source color/spin 

      psi = zero ; 

      StopWatch swatch;
      swatch.start();
      SystemSolverResults_t res = (*qprop)(psi, q_source);
      swatch.stop();
      double time_in_sec  = swatch.getTimeInSeconds();
      ncg_had += res.n_count;
      
      push(xml_out,"Qprop");
      write(xml_out, "Mass" , Mass);
      write(xml_out, "RsdCG", RsdCG);
      write(xml_out, "n_count", res.n_count);
      write(xml_out, "time_in_sec",time_in_sec );
      pop(xml_out); //Qprop
  
      FermToProp(psi, quark_propagator_4link, color_source);

    }


    staggered_pion_singlet pion_singlet(t_length,u);
    pion_singlet.compute(quark_propagator,quark_propagator_4link,j_decay);
    pion_singlet.dump(t_source,xml_out ) ;

    pop(xml_out);  //Computation_4link_pseudoscalar

    return ncg_had ;
  }
  /***************************************************************************/

  int compute_vary_singlet_ps(LatticeStaggeredFermion & psi,
	           LatticeStaggeredPropagator & quark_propagator_Lsink_Lsrc,
	       	   LatticeStaggeredPropagator & quark_propagator_Fsink_Lsrc,
		   LatticeStaggeredPropagator & quark_propagator_Lsink_Fsrc,
		   LatticeStaggeredPropagator & quark_propagator_Fsink_Fsrc,
		   stag_src_type type_of_src,
		   bool gauge_shift,
		   bool sym_shift,
		   const multi1d<LatticeColorMatrix> & u ,
		   const multi1d<LatticeColorMatrix> & u_smr ,
		   Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
		   XMLWriter & xml_out,
		   Real RsdCG, Real Mass, 
		   int j_decay, int t_source, int t_length, 
			      int fuzz_width){




    LatticeStaggeredFermion q_source_4link ;

    LatticeStaggeredFermion q_source ;
    LatticeStaggeredFermion q_source_fuzz ;
    int ncg_had = 0 ;

    LatticeStaggeredPropagator quark_propagator_4link ;

    q_source = zero ;

    // safety checks
    check_qprop_source_compatability(type_of_src, gauge_shift, sym_shift,
				     SINGLETS);

    if( gauge_shift ){

      if((type_of_src == LOCAL_SRC) || (type_of_src == FUZZED_SRC)){
	QDPIO::cerr << "Conflicting source and shift types in " <<
	  "util_compute_quark_prop_s.cc" <<endl;
	exit(0);
      }


    }else{

      if (sym_shift){
	QDPIO::cerr << "Conflicting source and shift types in " <<
	  "util_compute_singlet_s.cc" <<endl;
	QDPIO::cerr << "symmetric shifts not implemented for LLL " <<
	  "non-gauge-invariant sources" << endl;

	if(type_of_src == LOCAL_SRC){
	  QDPIO::cerr << "You asked for LOCAL_SRC, sym_shift" << endl;
	}
	if(type_of_src == FUZZED_SRC){
	  QDPIO::cerr << "You asked for FUZZED_SRC, sym_shift" << endl;
	}

	exit(0);
      }
    
    
      if(type_of_src == GAUGE_INVAR_LOCAL_SOURCE){
	QDPIO::cout << "Conflicting source and shift types in " <<
	  "util_compute_quark_prop_s.cc" <<endl;
	exit(0);
      }
    }


    /*** generate the 4-link quark propagator ****/
    push(xml_out,"Computation_4link_pseudoscalar");


    for(int color_source = 0; color_source < Nc; ++color_source) {

      if( type_of_src == LOCAL_SRC ){

	q_source = zero ;
	multi1d<int> coord(Nd);
	coord[0]=1; coord[1] = 1; coord[2] = 1; coord[3] = 1;
	srcfil(q_source, coord,color_source ) ;
      } else {
	if( type_of_src == GAUGE_INVAR_LOCAL_SOURCE  ) {

	  q_source = zero ;
	  multi1d<int> coord(Nd);
	  
	  // start with local source 
	  coord[0]= 0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
	  srcfil(q_source, coord,color_source ) ;


	  
	  // now do the shift
	  coord[0]=1; coord[1] = 1; coord[2] = 1; coord[3] = 1;
	
	  q_source_fuzz = q_source  ;
	  q_source = shiftDeltaPropCov(coord,q_source_fuzz,u,
				       sym_shift);

	}
      }

      // we fuzz the local quark propagator so it does not
      // make sense to fuzz the shited source quark propagator


      // Compute the propagator for given source color/spin 

      psi = zero ; 

      StopWatch swatch;
      swatch.start();
      SystemSolverResults_t res = (*qprop)(psi, q_source);
      swatch.stop();
      double time_in_sec  = swatch.getTimeInSeconds();
      ncg_had += res.n_count;
      
      push(xml_out,"Qprop");
      write(xml_out, "Mass" , Mass);
      write(xml_out, "RsdCG", RsdCG);
      write(xml_out, "n_count", res.n_count);
      write(xml_out, "time_in_sec",time_in_sec );
      pop(xml_out);
  
      FermToProp(psi, quark_propagator_4link, color_source);

    }


    staggered_pion_singlet pion_singlet_LsrcLsink(t_length,u);
    staggered_pion_singlet pion_singlet_FsrcLsink(t_length,u);
    staggered_pion_singlet pion_singlet_LsrcFsink(t_length,u);
    staggered_pion_singlet pion_singlet_FsrcFsink(t_length,u);
 
    push(xml_out,"LocalSource");
    push(xml_out,"LocalSink");
    pion_singlet_LsrcLsink.compute(quark_propagator_Lsink_Lsrc,
				   quark_propagator_4link,j_decay);
    pion_singlet_LsrcLsink.dump(t_source,xml_out ) ;
    pop(xml_out); //local sink


    push(xml_out,"FuzzSink");
    pion_singlet_LsrcFsink.compute(quark_propagator_Fsink_Lsrc,
				   quark_propagator_4link,j_decay);
    pion_singlet_LsrcFsink.dump(t_source,xml_out ) ;
    pop(xml_out); //fuzz sink
    pop(xml_out); //local source

    push(xml_out,"FuzzSource");
    push(xml_out,"LocalSink");
    pion_singlet_FsrcLsink.compute(quark_propagator_Lsink_Fsrc,
				   quark_propagator_4link,j_decay);
    pion_singlet_FsrcLsink.dump(t_source,xml_out ) ;

    pop(xml_out); //local sink

    push(xml_out,"FuzzSink");
    pion_singlet_FsrcFsink.compute(quark_propagator_Fsink_Fsrc,
				   quark_propagator_4link,j_decay);
    pion_singlet_FsrcFsink.dump(t_source,xml_out ) ;
    pop(xml_out); //fuzz sink

    pop(xml_out); //fuzz source


    pop(xml_out);  //Computation_4link_pseudoscalar

    return ncg_had ;



  }
} // end of namespace
