//
// Wrapper code to compute the connected part of the 
// pseudocalar singlet meson (4 link operator)
//
//

#include "handle.h"
#include "actions/ferm/fermbcs/fermbcs.h"
#include "actions/ferm/fermacts/fermacts_s.h"
#include "meas/hadron/hadron_s.h"
#include "meas/smear/fuzz_smear.h"
#include "meas/hadron/srcfil.h"
#include "util/ferm/transf.h"
#include "meas/hadron/pion_sing_s.h"


#include "util_compute_quark_prop_s.h"

namespace Chroma 
{ 


int compute_singlet_ps(LatticeStaggeredFermion & psi,
		       LatticeStaggeredPropagator quark_propagator,
		       stag_src_type type_of_src, 
		       const multi1d<LatticeColorMatrix> & u , 
		       Handle<const SystemSolver<LatticeStaggeredFermion> > & qprop,
		       XMLWriter & xml_out,
		       Real RsdCG, Real Mass, 
		       int j_decay, int t_source, int t_length)
{
  LatticeStaggeredFermion q_source ;
  LatticeStaggeredFermion q_source_fuzz ; 
  int ncg_had = 0 ;

  LatticeStaggeredPropagator quark_propagator_4link ;

  q_source = zero ;

  /*** generate the 4-link quark propagator ****/
  push(xml_out,"Computation_4link_pseudoscalar");

    for(int color_source = 0; color_source < Nc; ++color_source)
    {

      if( type_of_src == LOCAL_SRC )
	{
	  q_source = zero ;
	  multi1d<int> coord(Nd);
	  coord[0]=1; coord[1] = 1; coord[2] = 1; coord[3] = 1;
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
	  coord[0]=1; coord[1] = 1; coord[2] = 1; coord[3] = 1;
	  q_source_fuzz = q_source  ;
	  q_source = shiftDeltaPropCov(coord,q_source_fuzz,u,
				       false); 

	}

      // we fuzz the local quark propagator so it does not
      // make sense to fuzz the shited source quark propagator


      // Compute the propagator for given source color/spin 
      // int n_count;

      psi = zero ; 

      StopWatch swatch;
      swatch.start();
      int n_count = (*qprop)(psi, q_source);
      swatch.stop();
      double time_in_sec  = swatch.getTimeInSeconds();
      ncg_had += n_count;
      
      push(xml_out,"Qprop");
      write(xml_out, "Mass" , Mass);
      write(xml_out, "RsdCG", RsdCG);
      write(xml_out, "n_count", n_count);
      write(xml_out, "time_in_sec",time_in_sec );
      pop(xml_out);
  

      FermToProp(psi, quark_propagator_4link, color_source);

    }


    staggered_pion_singlet pion_singlet(t_length,u);
    pion_singlet.use_gauge_invar() ;
    pion_singlet.compute(quark_propagator,quark_propagator_4link,j_decay);
    pion_singlet.dump(t_source,xml_out ) ;

    pop(xml_out);

  return ncg_had ;
}



} // end of namespace
