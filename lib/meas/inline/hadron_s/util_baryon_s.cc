// $Id: util_baryon_s.cc,v 3.1 2009-05-08 11:09:21 mcneile Exp $
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
		       int j_decay, int tlength,
		       bool binary_baryon_dump, std::string filename)
{
  int bc_spec = 0 ;
  multi1d<int> coord(Nd);


  coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
  
  multi1d<Complex>  barprop(tlength) ;
  
  baryon_s(quark_propagator_a,quark_propagator_b,quark_propagator_c,
	   barprop,coord,j_decay, bc_spec) ;
  

  if(binary_baryon_dump){
      const int magic_number = 66618 ;
      BinaryFileWriter speedy ;
      speedy.open(filename);
      write(speedy,magic_number) ;
      write(speedy,tlength) ;
      write(speedy, barprop);;
  }else{
    write(xml_out, name, barprop);
  }
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


void compute_vary_baryon_s(XMLWriter &xml_out, int t_source, int fuzz_width,
			   int j_decay, int t_len, 
			   LatticeStaggeredPropagator & quark_propagator_Lsink_Lsrc,
			   LatticeStaggeredPropagator & quark_propagator_Fsink_Lsrc,
			   LatticeStaggeredPropagator & quark_propagator_Lsink_Fsrc,
			   LatticeStaggeredPropagator & quark_propagator_Fsink_Fsrc,
			   bool binary_baryon_dump, std::string binary_name)

{


  string filename_base;


      //
      // compute some simple baryon correlators
      //

      push(xml_out, "baryon_correlators");

      // describe the source
      string NN ;
      write(xml_out, "source_time", t_source);
      push(xml_out, "smearing_info");
      NN = "L" ; 
      write_smearing_info(NN, LOCAL_SRC,xml_out,fuzz_width) ;

      NN = "F" ; 
      write_smearing_info(NN,FUZZED_SRC,xml_out,fuzz_width) ;
    
      pop(xml_out);

      // write out the baryon correlators 
      string b_tag("srcLLL_sinkLLL_nucleon") ;  


      filename_base=binary_name+b_tag;
      ks_compute_baryon(b_tag,quark_propagator_Lsink_Lsrc, 
			quark_propagator_Lsink_Lsrc, 
			quark_propagator_Lsink_Lsrc, 
			xml_out, j_decay, 
			t_len,
			binary_baryon_dump,filename_base) ;

      // single quark fuzzed

      b_tag = "srcLLL_sinkFLL_nucleon" ;
      filename_base=binary_name+b_tag;
      ks_compute_baryon(b_tag,
			quark_propagator_Fsink_Lsrc, 
			quark_propagator_Lsink_Lsrc, 
			quark_propagator_Lsink_Lsrc, 
			xml_out, j_decay, 
			t_len,
			binary_baryon_dump,filename_base) ;

      b_tag = "srcFLL_sinkLLL_nucleon" ;
      filename_base=binary_name+b_tag;
      ks_compute_baryon(b_tag,
			quark_propagator_Lsink_Fsrc, 
			quark_propagator_Lsink_Lsrc, 
			quark_propagator_Lsink_Lsrc, 
			xml_out, j_decay, 
			t_len,
			binary_baryon_dump,filename_base) ;

      b_tag = "srcFLL_sinkFLL_nucleon" ;
      filename_base=binary_name+b_tag;
      ks_compute_baryon(b_tag,
			quark_propagator_Fsink_Fsrc, 
			quark_propagator_Lsink_Lsrc, 
			quark_propagator_Lsink_Lsrc, 
			xml_out, j_decay, 
			t_len,
			binary_baryon_dump,filename_base) ;



      // double quark fuzzed

      b_tag = "srcLLL_sinkFFL_nucleon" ;
      filename_base=binary_name+b_tag;
      ks_compute_baryon(b_tag,
			quark_propagator_Fsink_Lsrc, 
			quark_propagator_Fsink_Lsrc, 
			quark_propagator_Lsink_Lsrc, 
			xml_out, j_decay, 
			t_len,
			binary_baryon_dump,filename_base) ;

      b_tag = "srcFFL_sinkLLL_nucleon" ;
      filename_base=binary_name+b_tag;
      ks_compute_baryon(b_tag,
			quark_propagator_Lsink_Fsrc, 
			quark_propagator_Lsink_Fsrc, 
			quark_propagator_Lsink_Lsrc, 
			xml_out, j_decay, 
			t_len,
			binary_baryon_dump,filename_base) ;

      b_tag = "srcFFL_sinkFFL_nucleon" ;
      filename_base=binary_name+b_tag;
      ks_compute_baryon(b_tag,
			quark_propagator_Fsink_Fsrc, 
			quark_propagator_Fsink_Fsrc, 
			quark_propagator_Lsink_Lsrc, 
			xml_out, j_decay, 
			t_len, 
			binary_baryon_dump,filename_base) ;


      // treble quark fuzzed

      b_tag = "srcLLL_sinkFFF_nucleon" ;
      filename_base=binary_name+b_tag;
      ks_compute_baryon(b_tag,
			quark_propagator_Fsink_Lsrc, 
			quark_propagator_Fsink_Lsrc, 
			quark_propagator_Fsink_Lsrc, 
			xml_out, j_decay, 
			t_len,
			binary_baryon_dump,filename_base) ;

      b_tag = "srcFFF_sinkLLL_nucleon" ;
      filename_base=binary_name+b_tag;
      ks_compute_baryon(b_tag,
			quark_propagator_Lsink_Fsrc, 
			quark_propagator_Lsink_Fsrc, 
			quark_propagator_Lsink_Fsrc, 
			xml_out, j_decay, 
			t_len,
			binary_baryon_dump,filename_base) ;

      b_tag = "srcFFF_sinkFFF_nucleon" ;
      filename_base=binary_name+b_tag;
      ks_compute_baryon(b_tag,
			quark_propagator_Fsink_Fsrc, 
			quark_propagator_Fsink_Fsrc, 
			quark_propagator_Fsink_Fsrc, 
			xml_out, j_decay, 
			t_len,
			binary_baryon_dump,filename_base) ;



      pop(xml_out);  // baryon correlators



}




} // end of chroma namespace

