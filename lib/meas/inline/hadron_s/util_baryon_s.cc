// $Id: util_baryon_s.cc,v 1.1 2005-08-25 16:38:40 mcneile Exp $
/*! \file
 *  \brief Wrappers to compute staggered baryon correlators
 *
 */

#include "meas/hadron/baryon_s.h"
#include "util_compute_quark_prop_s.h"

//
// wrapper routine for baryon operators
//

namespace Chroma 
{ 


void ks_compute_baryon(string name,
		       LatticeStaggeredPropagator & quark_propagator, 
		       XMLWriter & xml_out, 
		       int j_decay, int tlength)
{
  int bc_spec = 0 ;
  multi1d<int> coord(Nd);
  coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
  
  multi1d<Complex>  barprop(tlength) ;
  
  baryon_s(quark_propagator,barprop,
	   coord,j_decay, bc_spec) ;
  

  write(xml_out, name, barprop);

 

}



//
// wrapper routine for baryon operators
//

void ks_compute_baryon(string name,
		       LatticeStaggeredPropagator & quark_propagator_a, 
		       LatticeStaggeredPropagator & quark_propagator_b, 
		       LatticeStaggeredPropagator & quark_propagator_c, 
		       XMLWriter & xml_out, 
		       int j_decay, int tlength)
{
  int bc_spec = 0 ;
  multi1d<int> coord(Nd);
  coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
  
  multi1d<Complex>  barprop(tlength) ;
  
  baryon_s(quark_propagator_a,quark_propagator_b,quark_propagator_c,
	   barprop,coord,j_decay, bc_spec) ;
  
  write(xml_out, name, barprop);

}



void write_smearing_info(string name, stag_src_type type_of_src,
			 XMLWriter &xml_out, int fuzz_width )
{


  push(xml_out, "smearing_basis");
  write(xml_out, "element_tag", name);

  if( type_of_src == LOCAL_SRC )
    { write(xml_out, "source_type", "LOCAL_SRC"); }
  else if( type_of_src == GAUGE_INVAR_LOCAL_SOURCE  )
    { write(xml_out, "source_type", "GAUGE_INVAR_LOCAL_SOURCE"); }
  else if( type_of_src == FUZZED_SRC )
    { 
      write(xml_out, "source_type", "FUZZED_SRC"); 
      write(xml_out, "fuzzed_width", fuzz_width); 
    }

  pop(xml_out);

}




} // end of chroma namespace

