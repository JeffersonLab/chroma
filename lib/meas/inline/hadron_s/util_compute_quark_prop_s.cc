// $Id: util_compute_quark_prop_s.cc,v 3.3 2008-06-26 15:20:10 mcneile Exp $
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

    }else{
      if( type_of_src == GAUGE_INVAR_LOCAL_SOURCE  ) {
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

      }else{

	QDPIO::cerr << "Conflicting source and shift types in " <<
	  "util_compute_quark_prop_s.cc" <<endl;
	QDPIO::cerr << "double-check your source --- no fuzz-smearing here" 
		    << endl;
	exit(0);

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





} // end of namespace
