//
// Wrapper code to compute a staggered quark propagtor.
//
//
//

#include "handle.h"
//#include "actions/ferm/fermbcs/fermbcs.h"
#include "actions/ferm/fermacts/fermacts_s.h"
#include "meas/hadron/hadron_s.h"
#include "meas/smear/fuzz_smear.h"
#include "meas/sources/srcfil.h"
#include "meas/sources/dilute_gauss_src_s.h"

#include "util_compute_quark_prop_s.h"

namespace Chroma { 

  enum func_flag{
    QPROP_FUZZ,
    QPROP_NO_FUZZ,
    SINGLETS
  };

  typedef func_flag func_flag_type;
  /************************************************************************/

  int check_qprop_source_compatability(stag_src_type type_of_src,
				       bool gauge_shift,
				       bool sym_shift, func_flag_type fflag){

    string calling_function;
    string str_source_type;
    string do_sym;
    string do_gauge;
    string more_info;

    bool problem = false;

    switch (fflag){
    case QPROP_FUZZ:
      calling_function = "compute_quark_propagator_s (fuzzed version)";
      break;
    case QPROP_NO_FUZZ:
      calling_function = "compute_quark_propagator_s (non-fuzzed version)";
      break;
    case SINGLETS:
      calling_function = "compute_singlet_ps";
      break;
    default:
      calling_function = "unknown function";
    }

    switch (type_of_src){
    case LOCAL_SRC:
      str_source_type="LOCAL_SRC";
      break;
    case GAUGE_INVAR_LOCAL_SOURCE:
      str_source_type="GAUGE_INVAR_LOCAL_SOURCE";
      break;
    case FUZZED_SRC:
      str_source_type="FUZZED_SRC";
      break;
    case NOISY_LOCAL_SOURCE:
      str_source_type="NOISY_LOCAL_SOURCE";
      break;
    case LOAD_IN_SOURCE:
      str_source_type="LOAD_IN_SOURCE";
      break;
    default:
      QDPIO::cerr << "Unknown type of source in" << str_source_type << endl;
      exit (1);
    }

    if(gauge_shift){
      do_gauge = "gauge_shift=true";
    }else{
      do_gauge = "gauge_shift=false";
    }

    if(sym_shift){
      do_sym = "sym_shift=true";
    }else{
      do_sym = "sym_shift=false";
    }


    if( gauge_shift ){

      //      if((type_of_src == LOCAL_SRC) || (type_of_src == FUZZED_SRC)){
      if((type_of_src == LOCAL_SRC)){
	problem = true;
	more_info = 
	  "Source type not gauge invariant.\nProbably not what you want.";
      }
    }else{
      if (sym_shift){
	problem = true;
	more_info = 
	  "Symmetric shifts not implemented for non-gauge-invariant sources";
      }
      if(type_of_src == GAUGE_INVAR_LOCAL_SOURCE){
	problem = true;
	more_info = "GAUGE_INVAR_LOCAL_SOURCE needs gauge shifting";
      }
    }

    if(fflag==QPROP_NO_FUZZ){
      if(type_of_src == FUZZED_SRC){
	problem = true;
	more_info = "Somehow called no-fuzzing version with FUZZED_SRC";
      }
    }

    if(problem){
      QDPIO::cerr << "In "<< calling_function <<" you used: "<< endl;
      QDPIO::cerr << "  "<< str_source_type<< endl;
      QDPIO::cerr << "  "<< do_gauge << endl;
      QDPIO::cerr << "  "<< do_sym << endl;
      QDPIO::cerr << more_info << endl;
    
      exit(1);
    }


    return 1;

  }


  /**************************************************************************/

  int compute_quark_propagator_s(LatticeStaggeredFermion & psi,
				 stag_src_type type_of_src,
				 bool gauge_shift,
				 bool sym_shift,
				 int fuzz_width,
				 const multi1d<LatticeColorMatrix> & u ,
				 multi1d<LatticeColorMatrix> & u_smr,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass, 
				 int j_decay,
				 int src_ind, int color_source){

    LatticeStaggeredFermion q_source ;
    LatticeStaggeredFermion q_source_fuzz ; 
    int ncg_had = 0 ;

static int flag=1;

    QDPIO::cout << "Inversion for Color =  " << color_source << endl;
    q_source = zero ;

    // safety checks
    check_qprop_source_compatability(type_of_src, gauge_shift, sym_shift,
				     QPROP_FUZZ);



    if( type_of_src == LOCAL_SRC ){
      q_source = zero ;
      multi1d<int> coord(Nd);

      PropIndexTodelta(src_ind, coord) ;
      srcfil(q_source, coord,color_source ) ;

    }else{
      if( type_of_src == GAUGE_INVAR_LOCAL_SOURCE  ) {
	q_source = zero ;
	multi1d<int> coord(Nd);

	// start with local source 
	coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;

	srcfil(q_source, coord,color_source ) ;

	// now do the shift
	PropIndexTodelta(src_ind, coord) ;
	q_source_fuzz = q_source  ;

	q_source = shiftDeltaPropCov(coord, q_source_fuzz, u, sym_shift);

	// sym-shift true for symmetric shifting!!!!!!!

      }else{
	if( type_of_src == FUZZED_SRC ){

	  q_source = zero ;
	  multi1d<int> coord(Nd);

	  PropIndexTodelta(src_ind, coord) ;
	  srcfil(q_source, coord,color_source ) ;

	  fuzz_smear(u_smr, q_source,q_source_fuzz, fuzz_width, j_decay) ;

	  q_source = q_source_fuzz  ;
	}
      }
    }


    // Use the last initial guess as the current guess

    // Compute the propagator for given source color/spin 
    // int n_count;

    StopWatch swatch;
    swatch.start();

    SystemSolverResults_t res = (*qprop)(psi, q_source);
    swatch.stop();
    double time_in_sec  = swatch.getTimeInSeconds();
    
    ncg_had += res.n_count;

    // this is done for xmldif reasons
    if( src_ind == 0 ){
      push(xml_out,"Qprop");
      write(xml_out, "Staggered_src_tag" , src_ind);
      write(xml_out, "Mass" , Mass);
      write(xml_out, "RsdCG", RsdCG);
      write(xml_out, "n_count", res.n_count);
      write(xml_out, "time_in_sec",time_in_sec );
      pop(xml_out);
    }

    return ncg_had ;
  }





  /**************************************************************************/

#if 0
  int compute_quark_propagator_s(LatticeStaggeredFermion & psi,
				 stag_src_type type_of_src,
				 bool gauge_shift,
				 bool sym_shift,
				 int fuzz_width,
				 const multi1d<LatticeColorMatrix> & u ,
				 multi1d<LatticeColorMatrix> & u_smr,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass, 
				 int j_decay,
				 int src_ind, int color_source, 
                                 LatticeStaggeredFermion & q_source_in ){

    LatticeStaggeredFermion q_source ;
    LatticeStaggeredFermion q_source_fuzz ; 
    int ncg_had = 0 ;

static int flag=1;

    QDPIO::cout << "Inversion for Color =  " << color_source << endl;
    q_source = zero ;

    // safety checks
    check_qprop_source_compatability(type_of_src, gauge_shift, sym_shift,
				     QPROP_FUZZ);



    if( type_of_src == LOCAL_SRC ){
      q_source = zero ;
      multi1d<int> coord(Nd);

      PropIndexTodelta(src_ind, coord) ;
      srcfil(q_source, coord,color_source ) ;

    }else{
      if( type_of_src == GAUGE_INVAR_LOCAL_SOURCE  ) {
	q_source = zero ;
	multi1d<int> coord(Nd);

	// start with local source 
	coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;

	srcfil(q_source, coord,color_source ) ;

	// now do the shift
	PropIndexTodelta(src_ind, coord) ;
	q_source_fuzz = q_source  ;

	q_source = shiftDeltaPropCov(coord, q_source_fuzz, u, sym_shift);

	// sym-shift true for symmetric shifting!!!!!!!

      }else{
	if( type_of_src == FUZZED_SRC ){

	  q_source = zero ;
	  multi1d<int> coord(Nd);

	  PropIndexTodelta(src_ind, coord) ;
	  srcfil(q_source, coord,color_source ) ;

	  fuzz_smear(u_smr, q_source,q_source_fuzz, fuzz_width, j_decay) ;

	  q_source = q_source_fuzz  ;
	}
      }
    }


	if( type_of_src == LOAD_IN_SOURCE  ){
	    q_source = q_source_in ; 
	}
	else
	{
	    q_source_in = q_source ; 
	}




    // Use the last initial guess as the current guess

    // Compute the propagator for given source color/spin 
    // int n_count;

    StopWatch swatch;
    swatch.start();

    SystemSolverResults_t res = (*qprop)(psi, q_source);
    swatch.stop();
    double time_in_sec  = swatch.getTimeInSeconds();
    
    ncg_had += res.n_count;

    // this is done for xmldif reasons
    if( src_ind == 0 ){
      push(xml_out,"Qprop");
      write(xml_out, "Staggered_src_tag" , src_ind);
      write(xml_out, "Mass" , Mass);
      write(xml_out, "RsdCG", RsdCG);
      write(xml_out, "n_count", res.n_count);
      write(xml_out, "time_in_sec",time_in_sec );
      pop(xml_out);
    }

    return ncg_had ;
  }

#endif

  /***************************************************************************/

  int compute_quark_propagator_s(LatticeStaggeredFermion & psi,
				 stag_src_type type_of_src,
				 bool gauge_shift,
				 bool sym_shift,
				 const multi1d<LatticeColorMatrix> & u ,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass, 
				 int j_decay,
				 int src_ind, int color_source, int t_source){

    LatticeStaggeredFermion q_source ;
    LatticeStaggeredFermion q_source_fuzz ; 
    int ncg_had = 0 ;

    QDPIO::cout << "Inversion for Color =  " << color_source << endl;
    q_source = zero ;

    // safety checks
    check_qprop_source_compatability(type_of_src, gauge_shift, sym_shift,
				     QPROP_NO_FUZZ);


    if( type_of_src == LOCAL_SRC ){
      q_source = zero ;
      multi1d<int> coord(Nd);

      PropIndexTodelta(src_ind, coord) ;
      srcfil(q_source, coord,color_source ) ;

    }
    else if( type_of_src == GAUGE_INVAR_LOCAL_SOURCE  ) {
      
	q_source = zero ;
	multi1d<int> coord(Nd);

	// start with local source 
	coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
        coord[j_decay] = t_source ;
	srcfil(q_source, coord,color_source ) ;

	// now do the shift
	PropIndexTodelta(src_ind, coord) ;
	q_source_fuzz = q_source  ;

	q_source = shiftDeltaPropCov(coord, q_source_fuzz, u, sym_shift);

	// sym-shift true for symmetric shifting!!!!!!!
    }
    else if( type_of_src ==  NOISY_LOCAL_SOURCE  ) {

      gaussian_color_src_on_slice(q_source, color_source,t_source, j_decay);  

      }
    else{

	QDPIO::cerr << "Conflicting source and shift types in " <<
	  "util_compute_quark_prop_s.cc" <<endl;
	QDPIO::cerr << "double-check your source --- no fuzz-smearing here" 
		    << endl;
	exit(0);

      }
  


    // Use the last initial guess as the current guess

    // Compute the propagator for given source color/spin 
    // int n_count;

    StopWatch swatch;
    swatch.start();

    SystemSolverResults_t res = (*qprop)(psi, q_source);
    swatch.stop();
    double time_in_sec  = swatch.getTimeInSeconds();
    
    ncg_had += res.n_count;

    // this is done for xmldif reasons
    if( src_ind == 0 ){
      push(xml_out,"Qprop");
      write(xml_out, "Staggered_src_tag" , src_ind);
      write(xml_out, "Mass" , Mass);
      write(xml_out, "RsdCG", RsdCG);
      write(xml_out, "n_count", res.n_count);
      write(xml_out, "time_in_sec",time_in_sec );
      pop(xml_out);
    }

    return ncg_had ;
  }




  int compute_quark_propagator_s(LatticeStaggeredFermion & psi,
				 stag_src_type type_of_src,
				 bool gauge_shift,
				 bool sym_shift,
				 const multi1d<LatticeColorMatrix> & u ,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass, 
				 int j_decay,
				 int src_ind, int color_source, int t_source,
				 LatticeStaggeredFermion & q_source_in ){

    LatticeStaggeredFermion q_source ;
    LatticeStaggeredFermion q_source_fuzz ; 
    int ncg_had = 0 ;

    QDPIO::cout << "Inversion for Color =  " << color_source << endl;
    q_source = zero ;

    // safety checks
    check_qprop_source_compatability(type_of_src, gauge_shift, sym_shift,
				     QPROP_NO_FUZZ);


    if( type_of_src == LOCAL_SRC ){
      q_source = zero ;
      multi1d<int> coord(Nd);

      PropIndexTodelta(src_ind, coord) ;
      srcfil(q_source, coord,color_source ) ;

    }
    else if( type_of_src == GAUGE_INVAR_LOCAL_SOURCE  ) {
      
	q_source = zero ;
	multi1d<int> coord(Nd);

	// start with local source 
	coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
        coord[j_decay] = t_source ;
	srcfil(q_source, coord,color_source ) ;

	// now do the shift
	PropIndexTodelta(src_ind, coord) ;
	q_source_fuzz = q_source  ;

	q_source = shiftDeltaPropCov(coord, q_source_fuzz, u, sym_shift);

	// sym-shift true for symmetric shifting!!!!!!!
    }
    else if( type_of_src ==  NOISY_LOCAL_SOURCE  ) {

      gaussian_color_src_on_slice(q_source, color_source,t_source, j_decay);  

      }
    else if( type_of_src == LOAD_IN_SOURCE  ){
	    q_source = q_source_in ; 
	}
    else{

	QDPIO::cerr << "Conflicting source and shift types in " <<
	  "util_compute_quark_prop_s.cc" <<endl;
	QDPIO::cerr << "double-check your source --- no fuzz-smearing here" 
		    << endl;
	cout << "DEBUG type_of_src = " << type_of_src << cout ;
	exit(0);

      }
  

    if( type_of_src != LOAD_IN_SOURCE  )
      {
	q_source_in = q_source ; 
      }



    // Use the last initial guess as the current guess

    // Compute the propagator for given source color/spin 
    // int n_count;

    StopWatch swatch;
    swatch.start();

    SystemSolverResults_t res = (*qprop)(psi, q_source);
    swatch.stop();
    double time_in_sec  = swatch.getTimeInSeconds();
    
    ncg_had += res.n_count;

    // this is done for xmldif reasons
    if( src_ind == 0 ){
      push(xml_out,"Qprop");
      write(xml_out, "Staggered_src_tag" , src_ind);
      write(xml_out, "Mass" , Mass);
      write(xml_out, "RsdCG", RsdCG);
      write(xml_out, "n_count", res.n_count);
      write(xml_out, "time_in_sec",time_in_sec );
      pop(xml_out);
    }

    return ncg_had ;
  }





  /**
    For non-generate inversion I need two quark propagators from the same
    source.

  **/


  int compute_quark_propagator_s(LatticeStaggeredFermion & psi1,
				 LatticeStaggeredFermion & psi2,
				 stag_src_type type_of_src,
				 bool gauge_shift,
				 bool sym_shift,
				 const multi1d<LatticeColorMatrix> & u ,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop1,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop2,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass1, Real Mass2, 
				 int j_decay,
				 int src_ind, int color_source, int t_source, 
                                 LatticeStaggeredFermion & q_source_in){

    LatticeStaggeredFermion q_source ;
    LatticeStaggeredFermion q_source_fuzz ; 
    int ncg_had = 0 ;
    string src_tag ;

    QDPIO::cout << "Inversion for Color =  " << color_source << endl;
    q_source = zero ;

    // safety checks
    check_qprop_source_compatability(type_of_src, gauge_shift, sym_shift,
				     QPROP_NO_FUZZ);


    if( type_of_src == LOCAL_SRC ){
      q_source = zero ;
      multi1d<int> coord(Nd);

      PropIndexTodelta(src_ind, coord) ;
      srcfil(q_source, coord,color_source ) ;
      src_tag = "LOCAL_SRC" ;
    }
    else if( type_of_src == GAUGE_INVAR_LOCAL_SOURCE  ) {
      
	q_source = zero ;
	multi1d<int> coord(Nd);

	// start with local source 
	coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
        coord[j_decay] = t_source ;
	srcfil(q_source, coord,color_source ) ;

	// now do the shift
	PropIndexTodelta(src_ind, coord) ;
	q_source_fuzz = q_source  ;

	q_source = shiftDeltaPropCov(coord, q_source_fuzz, u, sym_shift);

	// sym-shift true for symmetric shifting!!!!!!!
      src_tag = "GAUGE_INVAR_LOCAL_SOURCE" ;
    }
    else if( type_of_src ==  NOISY_LOCAL_SOURCE  ) {

      gaussian_color_src_on_slice(q_source, color_source,t_source, j_decay);  
      src_tag = "NOISY_LOCAL_SOURCE" ;
      }
    else if( type_of_src == LOAD_IN_SOURCE  ){
	    q_source = q_source_in ; 
	}

    else{
      src_tag = "SOURCE_ERROR" ;
	QDPIO::cerr << "Conflicting source and shift types in " <<
	  "util_compute_quark_prop_s.cc" <<endl;
	QDPIO::cerr << "double-check your source --- no fuzz-smearing here" 
		    << endl;
	cout << "DEBUG type_of_src = " << type_of_src << cout ;
	exit(0);

      }
  

    
    if( type_of_src != LOAD_IN_SOURCE  )
      {
	q_source_in = q_source ; 
      }


    // Use the last initial guess as the current guess

    // Compute the propagator for given source color/spin 
    // int n_count;

    StopWatch swatch;
    swatch.start();

    SystemSolverResults_t res1 = (*qprop1)(psi1, q_source);
    SystemSolverResults_t res2 = (*qprop2)(psi2, q_source);

    swatch.stop();
    double time_in_sec  = swatch.getTimeInSeconds();
    
    ncg_had += res1.n_count;
    ncg_had += res2.n_count;

    // this is done for xmldif reasons
    if( src_ind == 0 ){
      push(xml_out,"inversions");
      push(xml_out,"Qprop1");
      write(xml_out, "Staggered_src_tag" , src_tag);
      write(xml_out, "Mass" , Mass1);
      write(xml_out, "Final_RsdCG", res1.resid);
      write(xml_out, "n_count", res1.n_count);
      pop(xml_out);

      push(xml_out,"Qprop2");
      write(xml_out, "Staggered_src_tag" , src_tag);
      write(xml_out, "Mass" , Mass2);
      write(xml_out, "Final_RsdCG", res2.resid);
      write(xml_out, "n_count", res2.n_count);
      pop(xml_out);

      write(xml_out, "time_in_sec",time_in_sec );
      pop(xml_out);
    }

    return ncg_had ;
  }



  /**
    For non-degenerate inversion I need two quark propagators from the same
    source.

  **/


  int compute_quark_propagator_s(LatticeStaggeredFermion & psi1,
				 LatticeStaggeredFermion & psi2,
				 stag_src_type type_of_src,
				 bool gauge_shift,
				 bool sym_shift,
				 const multi1d<LatticeColorMatrix> & u ,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop1,
				 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop2,
				 XMLWriter & xml_out,
				 Real RsdCG, Real Mass1, Real Mass2, 
				 int j_decay,
				 int src_ind, int color_source, int t_source){

    LatticeStaggeredFermion q_source ;
    LatticeStaggeredFermion q_source_fuzz ; 
    int ncg_had = 0 ;
    string src_tag ;

    QDPIO::cout << "Inversion for Color =  " << color_source << endl;
    q_source = zero ;

    // safety checks
    check_qprop_source_compatability(type_of_src, gauge_shift, sym_shift,
				     QPROP_NO_FUZZ);


    if( type_of_src == LOCAL_SRC ){
      q_source = zero ;
      multi1d<int> coord(Nd);

      PropIndexTodelta(src_ind, coord) ;
      srcfil(q_source, coord,color_source ) ;
      src_tag = "LOCAL_SRC" ;
    }
    else if( type_of_src == GAUGE_INVAR_LOCAL_SOURCE  ) {
      
	q_source = zero ;
	multi1d<int> coord(Nd);

	// start with local source 
	coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
        coord[j_decay] = t_source ;
	srcfil(q_source, coord,color_source ) ;

	// now do the shift
	PropIndexTodelta(src_ind, coord) ;
	q_source_fuzz = q_source  ;

	q_source = shiftDeltaPropCov(coord, q_source_fuzz, u, sym_shift);

	// sym-shift true for symmetric shifting!!!!!!!
      src_tag = "GAUGE_INVAR_LOCAL_SOURCE" ;
    }
    else if( type_of_src ==  NOISY_LOCAL_SOURCE  ) {

      gaussian_color_src_on_slice(q_source, color_source,t_source, j_decay);  
      src_tag = "NOISY_LOCAL_SOURCE" ;
      }
    else{
      src_tag = "SOURCE_ERROR" ;
	QDPIO::cerr << "Conflicting source and shift types in " <<
	  "util_compute_quark_prop_s.cc" <<endl;
	QDPIO::cerr << "double-check your source --- no fuzz-smearing here" 
		    << endl;
	exit(0);

      }
  


    // Use the last initial guess as the current guess

    // Compute the propagator for given source color/spin 
    // int n_count;

    StopWatch swatch;
    swatch.start();

    SystemSolverResults_t res1 = (*qprop1)(psi1, q_source);
    SystemSolverResults_t res2 = (*qprop2)(psi2, q_source);

    swatch.stop();
    double time_in_sec  = swatch.getTimeInSeconds();
    
    ncg_had += res1.n_count;
    ncg_had += res2.n_count;

    // this is done for xmldif reasons
    if( src_ind == 0 ){
      push(xml_out,"inversions");
      push(xml_out,"Qprop1");
      write(xml_out, "Staggered_src_tag" , src_tag);
      write(xml_out, "Mass" , Mass1);
      write(xml_out, "Final_RsdCG", res1.resid);
      write(xml_out, "n_count", res1.n_count);
      pop(xml_out);

      push(xml_out,"Qprop2");
      write(xml_out, "Staggered_src_tag" , src_tag);
      write(xml_out, "Mass" , Mass2);
      write(xml_out, "Final_RsdCG", res2.resid);
      write(xml_out, "n_count", res2.n_count);
      pop(xml_out);

      write(xml_out, "time_in_sec",time_in_sec );
      pop(xml_out);
    }

    return ncg_had ;
  }





} // end of namespace
