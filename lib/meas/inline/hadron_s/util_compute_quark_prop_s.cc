//
// Wrapper code to compute a staggered quark propagtor.
//
//
//

#include "handle.h"
#include "actions/ferm/fermbcs/fermbcs.h"
#include "actions/ferm/fermacts/fermacts_s.h"
#include "meas/hadron/hadron_s.h"
#include "meas/smear/fuzz_smear.h"
#include "meas/sources/srcfil.h"

#include "util_compute_quark_prop_s.h"

namespace Chroma 
{ 


int compute_quark_propagator_s(LatticeStaggeredFermion & psi,
				stag_src_type type_of_src, 
				int fuzz_width,
				const multi1d<LatticeColorMatrix> & u , 
				multi1d<LatticeColorMatrix> & u_smr,
				Handle<const SystemSolver<LatticeStaggeredFermion> > & qprop,
				XMLWriter & xml_out,
				Real RsdCG, Real Mass, 
				int j_decay, 
				int src_ind, int color_source)
{
  LatticeStaggeredFermion q_source ;
  LatticeStaggeredFermion q_source_fuzz ; 
  int ncg_had = 0 ;

  QDPIO::cout << "Inversion for Color =  " << color_source << endl;
  q_source = zero ;

  if( type_of_src == LOCAL_SRC )
    {
      q_source = zero ;
      multi1d<int> coord(Nd);

      PropIndexTodelta(src_ind, coord) ; 
      srcfil(q_source, coord,color_source ) ;
    }
	  else if( type_of_src == GAUGE_INVAR_LOCAL_SOURCE  )
	    {
	      q_source = zero ;
	      multi1d<int> coord(Nd);
	      
	      // start with local source 
	      coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
	      srcfil(q_source, coord,color_source ) ;
	      
	      // now do the shift
	      PropIndexTodelta(src_ind, coord) ; 
	      q_source_fuzz = q_source  ;
	      q_source = shiftDeltaPropCov(coord,q_source_fuzz,u,
					   false); 

	    }
	  else if( type_of_src == FUZZED_SRC )
	    {
	      q_source = zero ;
	      multi1d<int> coord(Nd);

	      PropIndexTodelta(src_ind, coord) ; 
	      srcfil(q_source, coord,color_source ) ;


	      fuzz_smear(u_smr, q_source,q_source_fuzz, 
			 fuzz_width, j_decay) ; 

	      q_source = q_source_fuzz  ;
	    }



	  // Use the last initial guess as the current guess

	  // Compute the propagator for given source color/spin 
	  // int n_count;

        StopWatch swatch;
	swatch.start();

	int n_count = (*qprop)(psi, q_source);
	  swatch.stop();
	  double time_in_sec  = swatch.getTimeInSeconds();

    
	  ncg_had += n_count;

	// this is done for xmldif reasons
	  if( src_ind == 0 )
	  {
	    push(xml_out,"Qprop");
	    write(xml_out, "Staggered_src_tag" , src_ind);
	    write(xml_out, "Mass" , Mass);
	    write(xml_out, "RsdCG", RsdCG);
	    write(xml_out, "n_count", n_count);
	    write(xml_out, "time_in_sec",time_in_sec );
	    pop(xml_out);
	  }


	return ncg_had ;
}



} // end of namespace
