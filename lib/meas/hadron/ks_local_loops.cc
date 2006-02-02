/* + */
/* $Id: ks_local_loops.cc,v 2.2 2006-02-02 16:42:15 egregory Exp $ ($Date: 2006-02-02 16:42:15 $) */


#include "fermact.h"
#include "meas/hadron/ks_local_loops.h"
#include "meas/hadron/hadron_s.h"
#include "meas/sources/z2_src.h"

#include "meas/smear/fuzz_smear.h"
#include "meas/sources/dilute_gauss_src_s.h"

namespace Chroma {



void ks_local_loops(
		 Handle<const SystemSolver<LatticeStaggeredFermion> > & qprop,
		 LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi ,
		 const multi1d<LatticeColorMatrix> & u,
		 XMLWriter & xml_out, 
		 bool gauge_shift,
		 bool sym_shift,
		 int t_length,
		 Real Mass,
		 int Nsamp,
		 Real RsdCG,
		 int CFGNO,
		 int volume_source,
		 int src_seperation,
		 int j_decay){


    push(xml_out,"local_loops_s");

    // write common parameters
    write(xml_out, "Mass" , Mass);
    if( volume_source == Z2NOISE  ){
      write(xml_out, "Random_volume_source" , "Z2NOISE");
    }else if( volume_source == GAUSSIAN ){
      write(xml_out, "Random_volume_source" , "GAUSSIAN");

    }










    write(xml_out, "Number_of_samples" , Nsamp);


    ///////
    Stag_shift_option type_of_shift ; 
    if( gauge_shift ){
      if(sym_shift){
	type_of_shift = SYM_GAUGE_INVAR ;
      }else{
	type_of_shift = GAUGE_INVAR ;
      }
    }else{
      if(sym_shift){
	type_of_shift = SYM_NON_GAUGE_INVAR ;
      }else{
	type_of_shift = NON_GAUGE_INVAR ;
      }
    }
    //////

    // set up the loop code
    local_scalar_loop scalar_one_loop(t_length,Nsamp,
				      u,type_of_shift) ; 
    non_local_scalar_loop scalar_two_loop(t_length,Nsamp,
					  u,type_of_shift) ; 
    threelink_pseudoscalar_loop eta3_loop(t_length,Nsamp,
					  u,type_of_shift) ; 
    fourlink_pseudoscalar_loop eta4_loop(t_length,Nsamp,
					 u,type_of_shift) ; 


    // Seed the RNG with the cfg number for now
    QDP::Seed seed;
    seed = CFGNO;
    RNG::setrn(seed);

  int src_tslice=0;
  int src_color_ind = 0;
  int src_parity_ind = 0;
  int src_corner_ind =0;

  double coverage_fraction;

  for(int i = 0; i < Nsamp; ++i){
    psi = zero;   // note this is ``zero'' and not 0
    RNG::savern(seed);

    // Fill the volume with random noise 
    if( volume_source == GAUSSIAN  ){
      gaussian(q_source);
    }else if( volume_source == Z2NOISE ){
       z2_src(q_source); 
    }else if( volume_source == T_DILUTE_GAUSS ){
      printf("tslice=%d\n",src_tslice);
      gaussian_on_timeslice(q_source,src_tslice,j_decay);
      src_tslice++;
      if(src_tslice>=t_length){
	src_tslice=0;
      }
      coverage_fraction=1.0/t_length;
    }else if( volume_source == C_DILUTE_GAUSS ){ 
      printf("src_color=%d\n",src_color_ind);
      gaussian_color_src(q_source,src_color_ind);
      src_color_ind++;
      if(src_color_ind>=Nc){
	src_color_ind=0;
      }
      coverage_fraction=1.0/Nc;
    }else if( volume_source == P_DILUTE_GAUSS ){ 
      printf("src_parity=%d\n",src_parity_ind);
      gaussian_on_parity(q_source,src_parity_ind);
      src_parity_ind++;
      if(src_parity_ind > 1){
	src_parity_ind=0;
      }
      coverage_fraction=1.0/2;
    }else if( volume_source == CT_DILUTE_GAUSS ){
      printf("tslice=%d  src_color=%d\n",src_tslice,src_color_ind);
      gaussian_color_src_on_slice(q_source,src_color_ind,
				  src_tslice, j_decay);
      src_tslice++;
      if(src_tslice>=t_length){
	src_tslice=0;
	src_color_ind++;
	if(src_color_ind>=Nc){
	  src_color_ind=0;
	}
      }
      coverage_fraction=1.0/(Nc*t_length);
    }else if( volume_source == CP_DILUTE_GAUSS ){
      printf("src_color=%d  parity=%d  \n",src_color_ind, src_parity_ind);
      gaussian_color_src_on_parity(q_source, src_color_ind, src_parity_ind);
      src_parity_ind++;
      if(src_parity_ind > 1){
	src_parity_ind=0;
	src_color_ind++;
	if(src_color_ind>=Nc){
	  src_color_ind=0;
	}
      }
      coverage_fraction=1.0/(2*Nc);
    }else if( volume_source == PT_DILUTE_GAUSS ){
      printf("src_tslice=%d  parity=%d  \n",src_tslice, src_parity_ind);
      gaussian_parity_src_on_slice(q_source, src_parity_ind, src_tslice,
				   j_decay);
      src_parity_ind++;
      if(src_parity_ind > 1){
	src_parity_ind=0;
	src_tslice++;
	if(src_tslice>=t_length){
	  src_tslice=0;
	}
      }
      coverage_fraction=1.0/(2*t_length);
    }else if( volume_source == MOD_T_DILUTE_GAUSS ){
      printf("src_tslice=%d seperation=%d\n",src_tslice, src_seperation);
      gaussian_on_mod_timeslice(q_source, src_tslice, j_decay,
				src_seperation);
      src_tslice++;
      if(src_tslice>=t_length){
	src_tslice=0;
      }
      if(t_length%src_seperation==0){
	coverage_fraction=(1.0/(src_seperation));
      }else{
	coverage_fraction=9999999;
      }
    }else if( volume_source == CORNER_DILUTE_GAUSS ){
      printf("src_corner=%d \n",src_corner_ind);
      gaussian_on_corner(q_source, src_corner_ind);
      src_corner_ind++;
      if(src_corner_ind>=16){
	src_corner_ind=0;
      }
      coverage_fraction=1.0/16;
    }else if( volume_source == COR_DBL_T_DILUTE_GAUSS ){
      //DEBUG
       printf("i=%d src_corner_ind=%d src_tslice=%d\n",
      	     i,src_corner_ind,src_tslice);
      //fflush(stdout);
      //DEBUG

      gaussian_corner_on_dbl_slice(q_source, src_corner_ind, src_tslice,
				   j_decay);
      src_corner_ind++;
      if(src_corner_ind>=16){
	src_corner_ind=0;
	src_tslice++;
	if(src_tslice>=t_length){
	  src_tslice=0;
	}
      }
      coverage_fraction=1.0/(16*t_length);
    }else if( volume_source ==   COR_MOD_DBL_T_DILUTE_GAUSS ){
      //DEBUG
      printf("i=%d src_corner_ind=%d src_tslice=%d (cmdtdg)\n",
	     i,src_corner_ind,src_tslice);
      fflush(stdout);
      //DEBUG

      gaussian_corner_on_mod_dbl_slice(q_source,src_corner_ind, src_tslice,
				       j_decay, src_seperation);
      src_corner_ind++;
      if(src_corner_ind>=16){
	src_corner_ind=0;
	src_tslice++;
	if(src_tslice>=t_length){
	  src_tslice=0;
	}
      }
      if(t_length%src_seperation==0){
	coverage_fraction=1.0/(2*Nc);
      }else{
	coverage_fraction=9999999;
      }
    }else if( volume_source == C_MOD_T_DILUTE_GAUSS ){

      printf("tslice=%d  color=%d\n",src_tslice,src_color_ind);

      gaussian_color_src_on_mod_slice(q_source,src_color_ind,
				  src_tslice, j_decay, src_seperation);
      src_tslice++;
      if(src_tslice>=t_length){
	src_tslice=0;
	src_color_ind++;
	if(src_color_ind>=Nc){
	  src_color_ind=0;
	}
      }

      coverage_fraction=1.0/(Nc*t_length);
    }else{
      QDP_error_exit("Wrong type of volume source");
    }






    // Compute the solution vector for the particular source
    int n_count = (*qprop)(psi, q_source);
      
    push(xml_out,"Qprop_noise");
    write(xml_out, "Noise_number" , i);
    write(xml_out, "RsdCG" , RsdCG);
    write(xml_out, "n_count", n_count);
    write(xml_out, "Seed" , seed);
    pop(xml_out);


    scalar_one_loop.compute(q_source,psi,i) ;
    scalar_two_loop.compute(q_source,psi,i) ;
    eta3_loop.compute(q_source,psi,i) ;
    eta4_loop.compute(q_source,psi,i) ;

  } // Nsamples


  // write output from the 
  scalar_one_loop.dump(xml_out) ;
  scalar_two_loop.dump(xml_out) ;
  eta3_loop.dump(xml_out) ;
  eta4_loop.dump(xml_out) ;

  // end of this section
  pop(xml_out);

}




  //
  //  version used in test code.
  //

void ks_local_loops(
		 Handle<const SystemSolver<LatticeStaggeredFermion> > & qprop,
		 LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi ,
		 multi1d<LatticeColorMatrix> & u,
		 XMLFileWriter & xml_out, 
		 XMLReader & xml_in ,
		 int t_length,
		 Real Mass,
		 int Nsamp,
		 Real RsdCG,
		 int CFGNO,
		 int volume_source
		 )
{

    push(xml_out,"ks_local_loops");

    // write common parameters
    write(xml_out, "Mass" , Mass);
    if( volume_source == Z2NOISE  )
      write(xml_out, "Random_volume_source" , "Z2NOISE");
    else if( volume_source == GAUSSIAN )
      write(xml_out, "Random_volume_source" , "GAUSSIAN");

    write(xml_out, "Number_of_samples" , Nsamp);

  //
  //  parse input files
  //

  // the wrapped disconnected loops
    bool gauge_shift ;
    bool sym_shift ;
    try{
      read(xml_in, "/propagator/param/use_gauge_invar_oper", gauge_shift ) ;
    }catch (const string& e){
      QDPIO::cerr << "Error reading data: " << e << endl;
      throw;
    }
    try{
      read(xml_in, "/propagator/param/use_sym_shift_oper", sym_shift ) ;
    }catch (const string& e){
      QDPIO::cerr << "Error reading data: " << e << endl;
      throw;
    }


     ///////
    Stag_shift_option type_of_shift ; 
    if( gauge_shift ){
      if(sym_shift){
	type_of_shift = SYM_GAUGE_INVAR ;
      }else{
	type_of_shift = GAUGE_INVAR ;
      }
    }else{
      if(sym_shift){
	type_of_shift = SYM_NON_GAUGE_INVAR ;
      }else{
	type_of_shift = NON_GAUGE_INVAR ;
      }
    }
    //////


    // set up the loop code
    local_scalar_loop scalar_one_loop(t_length,Nsamp,
				      u,type_of_shift) ; 
    non_local_scalar_loop scalar_two_loop(t_length,Nsamp,
					  u,type_of_shift) ; 
    threelink_pseudoscalar_loop eta3_loop(t_length,Nsamp,
					  u,type_of_shift) ; 
    fourlink_pseudoscalar_loop eta4_loop(t_length,Nsamp,
					 u,type_of_shift) ; 


    // Seed the RNG with the cfg number for now
    QDP::Seed seed;
    seed = CFGNO;
    RNG::setrn(seed);


  for(int i = 0; i < Nsamp; ++i){
    psi = zero;   // note this is ``zero'' and not 0
    RNG::savern(seed);

    // Fill the volume with random noise 
    if( volume_source == GAUSSIAN  )
      gaussian(q_source);
    else if( volume_source == Z2NOISE )
      { z2_src(q_source); }

    // Compute the solution vector for the particular source
    int n_count = (*qprop)(psi, q_source);
      
    push(xml_out,"Qprop_noise");
    write(xml_out, "Noise_number" , i);
    write(xml_out, "RsdCG" , RsdCG);
    write(xml_out, "n_count", n_count);
    write(xml_out, "Seed" , seed);
    pop(xml_out);


    scalar_one_loop.compute(q_source,psi,i) ;
    scalar_two_loop.compute(q_source,psi,i) ;
    eta3_loop.compute(q_source,psi,i) ;
    eta4_loop.compute(q_source,psi,i) ;

  } // Nsamples


  // write output from the 
  scalar_one_loop.dump(xml_out) ;
  scalar_two_loop.dump(xml_out) ;
  eta3_loop.dump(xml_out) ;
  eta4_loop.dump(xml_out) ;

  // end of this section
  pop(xml_out);

}



//
//  fuzz the loops
//
//




void ks_fuzz_loops(
		 Handle<const SystemSolver<LatticeStaggeredFermion> > & qprop,
		 LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi ,
		 LatticeStaggeredFermion & psi_fuzz ,
		 const multi1d<LatticeColorMatrix> & u,
		 const multi1d<LatticeColorMatrix> & u_smr,
		 XMLWriter & xml_out, 
		 bool gauge_shift,
		 bool sym_shift,
		 int t_length,
		 Real Mass,
		 int Nsamp,
		 Real RsdCG,
		 int CFGNO,
		 int volume_source,
		 int fuzz_width, 
		 int j_decay
		 )
{


    push(xml_out,"fuzz_loops_s");

    // write common parameters
    write(xml_out, "Mass" , Mass);
    if( volume_source == Z2NOISE  )
      write(xml_out, "Random_volume_source" , "Z2NOISE");
    else if( volume_source == GAUSSIAN )
      write(xml_out, "Random_volume_source" , "GAUSSIAN");

    write(xml_out, "Number_of_samples" , Nsamp);
    write(xml_out, "fuzz_width" , fuzz_width);

      ///////
    Stag_shift_option type_of_shift ; 
    if( gauge_shift ){
      if(sym_shift){
	type_of_shift = SYM_GAUGE_INVAR ;
      }else{
	type_of_shift = GAUGE_INVAR ;
      }
    }else{
      if(sym_shift){
	type_of_shift = SYM_NON_GAUGE_INVAR ;
      }else{
	type_of_shift = NON_GAUGE_INVAR ;
      }
    }
    //////


   // set up the loop code
    local_scalar_loop scalar_one_loop(t_length,Nsamp,
				      u,type_of_shift) ; 
    non_local_scalar_loop scalar_two_loop(t_length,Nsamp,
					  u,type_of_shift) ; 
    threelink_pseudoscalar_loop eta3_loop(t_length,Nsamp,
					  u,type_of_shift) ; 
    fourlink_pseudoscalar_loop eta4_loop(t_length,Nsamp,
					 u,type_of_shift) ; 


    // Seed the RNG with the cfg number for now
    QDP::Seed seed;
    seed = CFGNO;
    RNG::setrn(seed);


  for(int i = 0; i < Nsamp; ++i){
    psi = zero;   // note this is ``zero'' and not 0
    RNG::savern(seed);

    // Fill the volume with random noise 
    if( volume_source == GAUSSIAN  )
      gaussian(q_source);
    else if( volume_source == Z2NOISE )
      { z2_src(q_source); }

    // Compute the solution vector for the particular source
    int n_count = (*qprop)(psi, q_source);
      
    push(xml_out,"Qprop_noise");
    write(xml_out, "Noise_number" , i);
    write(xml_out, "RsdCG" , RsdCG);
    write(xml_out, "n_count", n_count);
    write(xml_out, "Seed" , seed);
    pop(xml_out);


    fuzz_smear(u_smr, psi,psi_fuzz,
	       fuzz_width, j_decay) ;


    scalar_one_loop.compute(q_source,psi_fuzz,i) ;
    scalar_two_loop.compute(q_source,psi_fuzz,i) ;
    eta3_loop.compute(q_source,psi_fuzz,i) ;
    eta4_loop.compute(q_source,psi_fuzz,i) ;

  } // Nsamples


  // write output from the 
  scalar_one_loop.dump(xml_out) ;
  scalar_two_loop.dump(xml_out) ;
  eta3_loop.dump(xml_out) ;
  eta4_loop.dump(xml_out) ;

  // end of this section
  pop(xml_out);

}





}  // end namespace Chroma
