/* + */
/* $Id: ks_local_loops.cc,v 3.5 2008/02/24 11:29:36 mcneile Exp $ ($Date: 2008/02/24 11:29:36 $) */


#include "fermact.h"
#include "meas/hadron/ks_local_loops.h"
#include "meas/hadron/hadron_s.h"
#include "meas/sources/z2_src.h"
#include "meas/sources/srcfil.h"

#include "meas/smear/fuzz_smear.h"
#include "meas/sources/dilute_gauss_src_s.h"

namespace Chroma {



  void write_out_source_type(XMLWriter & xml_out, VolSrc_type volume_source){


    if( volume_source == Z2NOISE  ){
      write(xml_out, "Random_volume_source" , "Z2NOISE");

    }else if( volume_source == GAUSSIAN ){
      write(xml_out, "Random_volume_source" , "GAUSSIAN");

    }else if( volume_source == T_DILUTE_GAUSS ){
      write(xml_out, "Random_volume_source" , "T_DILUTE_GAUSS");

    }else if( volume_source == C_DILUTE_GAUSS ){ 
      write(xml_out, "Random_volume_source" , "C_DILUTE_GAUSS");

    }else if( volume_source == P_DILUTE_GAUSS ){ 
      write(xml_out, "Random_volume_source" , "P_DILUTE_GAUSS");

    }else if( volume_source == CT_DILUTE_GAUSS ){
      write(xml_out, "Random_volume_source" , "CT_DILUTE_GAUSS");

    }else if( volume_source == CP_DILUTE_GAUSS ){
      write(xml_out, "Random_volume_source" , "CP_DILUTE_GAUSS");

    }else if( volume_source == PT_DILUTE_GAUSS ){
      write(xml_out, "Random_volume_source" , "PT_DILUTE_GAUSS");

    }else if( volume_source == MOD_T_DILUTE_GAUSS ){
     write(xml_out, "Random_volume_source" , "MOD_T_DILUTE_GAUSS");

    }else if( volume_source == CORNER_DILUTE_GAUSS ){
     write(xml_out, "Random_volume_source" , "CORNER_DILUTE_GAUSS");

    }else if( volume_source == COR_DBL_T_DILUTE_GAUSS ){
     write(xml_out, "Random_volume_source" , "COR_DBL_T_DILUTE_GAUSS");

    }else if( volume_source ==   COR_MOD_DBL_T_DILUTE_GAUSS ){
     write(xml_out, "Random_volume_source" , "COR_MOD_DBL_T_DILUTE_GAUSS");

    }else if( volume_source == C_MOD_T_DILUTE_GAUSS ){
     write(xml_out, "Random_volume_source" , "C_MOD_T_DILUTE_GAUSS");
    }


  }


  /**************************************************************************/

  Real fill_volume_source(LatticeStaggeredFermion & q_source, 
			  VolSrc_type volume_source, int t_length, 
			  int *p_src_tslice, int *p_src_color_ind, 
			  int *p_src_parity_ind, int *p_src_corner_ind, 
			  int src_seperation, int j_decay){
 
  // Fill the volume with random noise 
 
    Real coverage_fraction;

   if( volume_source == GAUSSIAN  ){
      gaussian(q_source);
    }else if( volume_source == Z2NOISE ){
       z2_src(q_source); 

    }else if( volume_source == T_DILUTE_GAUSS ){
      gaussian_on_timeslice(q_source,(*p_src_tslice),j_decay);
      (*p_src_tslice)++;
      if((*p_src_tslice)>=t_length){
	(*p_src_tslice)=0;
      }
      coverage_fraction=1.0/t_length;

    }else if( volume_source == C_DILUTE_GAUSS ){ 
      gaussian_color_src(q_source,(*p_src_color_ind));
      (*p_src_color_ind)++;
      if((*p_src_color_ind) >= Nc){
	(*p_src_color_ind)=0;
      }
      coverage_fraction=1.0/Nc;

    }else if( volume_source == P_DILUTE_GAUSS ){ 
      gaussian_on_parity(q_source,(*p_src_parity_ind));
      (*p_src_parity_ind)++;
      if((*p_src_parity_ind) > 1){
	(*p_src_parity_ind) = 0;
      }
      coverage_fraction=1.0/2;

    }else if( volume_source == CT_DILUTE_GAUSS ){
      gaussian_color_src_on_slice(q_source,(*p_src_color_ind),
				  (*p_src_tslice), j_decay);
      (*p_src_tslice)++;
      if((*p_src_tslice)>=t_length){
	(*p_src_tslice)=0;
	(*p_src_color_ind)++;
	if((*p_src_color_ind) >=Nc){
	  (*p_src_color_ind)=0;
	}
      }
      coverage_fraction=1.0/(Nc*t_length);

    }else if( volume_source == CP_DILUTE_GAUSS ){
      gaussian_color_src_on_parity(q_source, (*p_src_color_ind), 
				   (*p_src_parity_ind));
      (*p_src_parity_ind)++;
      if((*p_src_parity_ind) > 1){
	(*p_src_parity_ind)=0;
	(*p_src_color_ind)++;
	if((*p_src_color_ind)>=Nc){
	  (*p_src_color_ind)=0;
	}
      }
      coverage_fraction=1.0/(2*Nc);

    }else if( volume_source == PT_DILUTE_GAUSS ){
      gaussian_parity_src_on_slice(q_source, (*p_src_parity_ind), 
				   (*p_src_tslice), j_decay);
      (*p_src_parity_ind)++;
      if((*p_src_parity_ind) > 1){
	(*p_src_parity_ind)=0;
	(*p_src_tslice)++;
	if((*p_src_tslice) >= t_length){
	  (*p_src_tslice) = 0;
	}
      }
      coverage_fraction=1.0/(2*t_length);

    }else if( volume_source == MOD_T_DILUTE_GAUSS ){
      gaussian_on_mod_timeslice(q_source, (*p_src_tslice), j_decay,
				src_seperation);
      (*p_src_tslice)++;
      if((*p_src_tslice)>=t_length){
	(*p_src_tslice)=0;
      }
      if(t_length%src_seperation==0){
	coverage_fraction=(1.0/(src_seperation));
      }else{
	coverage_fraction=9999999;
      }

    }else if( volume_source == CORNER_DILUTE_GAUSS ){
      gaussian_on_corner(q_source, (*p_src_corner_ind));
      (*p_src_corner_ind)++;
      if((*p_src_corner_ind)>=16){
	(*p_src_corner_ind)=0;
      }
      coverage_fraction=1.0/16;

    }else if( volume_source == COR_DBL_T_DILUTE_GAUSS ){
       gaussian_corner_on_dbl_slice(q_source, (*p_src_corner_ind), 
				   (*p_src_tslice), j_decay);
      (*p_src_corner_ind)++;
      if((*p_src_corner_ind) >= 16){
	(*p_src_corner_ind) = 0;
	(*p_src_tslice)++;
	if((*p_src_tslice)>=t_length){
	  (*p_src_tslice) = 0;
	}
      }
      coverage_fraction=1.0/(16*t_length);
 
    }else if( volume_source ==   COR_MOD_DBL_T_DILUTE_GAUSS ){
       gaussian_corner_on_mod_dbl_slice(q_source,(*p_src_corner_ind), 
					(*p_src_tslice), j_decay, 
					src_seperation);
      (*p_src_corner_ind)++;
      if((*p_src_corner_ind) >= 16){
	(*p_src_corner_ind)=0;
	(*p_src_tslice)++;
	if((*p_src_tslice)>=t_length){
	  (*p_src_tslice)=0;
	}
      }
      if(t_length%src_seperation==0){
	coverage_fraction=1.0/(2*Nc);
      }else{
	coverage_fraction=9999999;
      }

    }else if( volume_source == C_MOD_T_DILUTE_GAUSS ){
      gaussian_color_src_on_mod_slice(q_source,(*p_src_color_ind),
				  (*p_src_tslice), j_decay, src_seperation);
      (*p_src_tslice)++;
      if((*p_src_tslice)>=t_length){
	(*p_src_tslice)=0;
	(*p_src_color_ind)++;
	if((*p_src_color_ind)>=Nc){
	  (*p_src_color_ind) = 0;
	}
      }

      coverage_fraction=1.0/(Nc*t_length);

    }else{
      QDP_error_exit("Wrong type of volume source");
    }

   /**********************************************/
   /* local src for testing */

   //    q_source=zero;
   //   multi1d<int> coord(Nd);
   //   coord[0]=0; coord[1] = 0; coord[2] = 0; coord[3] = 0;
   //   srcfil(q_source, coord,0 ) ;
   //   srcfil(q_source, coord,1 ) ;
   //   srcfil(q_source, coord,2 ) ;





   /***********************************************/



   return coverage_fraction;


  }


  /**************************************************************************/



void ks_local_loops(
		 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
		 LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi ,
		 const multi1d<LatticeColorMatrix> & u,
		 XMLWriter & xml_out, 
		 bool gauge_shift,
		 bool sym_shift,
		 bool loop_checkpoint,
		 int t_length,
		 Real Mass,
		 int Nsamp,
		 Real RsdCG,
		 int CFGNO,
		 VolSrc_type volume_source,
		 int src_seperation,
		 int j_decay){


    push(xml_out,"local_loops_s");

   // write common parameters
    write(xml_out, "Mass" , Mass);

    write_out_source_type(xml_out, volume_source);

    write(xml_out, "Number_of_samples" , Nsamp);


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


    // set up the loop code
    local_scalar_loop                 scalar_one_loop(t_length,Nsamp,
						      u,type_of_shift) ;
    local_scalar_kilcup_loop          scalar_one_kilcup_loop(t_length, 
							     Nsamp, u,
							     type_of_shift);
    non_local_scalar_loop             scalar_two_loop(t_length,Nsamp,
						      u,type_of_shift) ;
    fourlink_scalar_loop scalar_four_loop(t_length, Nsamp,
					  u, type_of_shift);
    fourlink_scalar_kilcup_loop scalar_four_kilcup_loop(t_length, Nsamp,
							u, type_of_shift);
    threelink_pseudoscalar_loop       eta3_loop(t_length,Nsamp,
						u,type_of_shift) ;
    fourlink_pseudoscalar_loop        eta4_loop(t_length,Nsamp,
						u,type_of_shift) ;

    fourlink_pseudoscalar_kilcup_loop eta4_kilcup_loop(t_length,Nsamp,
						u,type_of_shift) ;

    zerolink_pseudoscalar_loop        eta0_loop(t_length, Nsamp, 
						u,type_of_shift) ;


    // Seed the RNG with the cfg number for now
    QDP::Seed seed;
    seed = CFGNO;
    RNG::setrn(seed);

  int src_tslice=0;
  int src_color_ind = 0;
  int src_parity_ind = 0;
  int src_corner_ind =0;

  Real coverage_fraction;

  for(int i = 0; i < Nsamp; ++i){
    psi = zero;   // note this is ``zero'' and not 0
    RNG::savern(seed);
    QDPIO::cout << "SEED = " << seed << endl;

    QDPIO::cout << "Noise sample: " << i << endl;

    // Fill the volume with random noise 
    coverage_fraction = fill_volume_source(q_source, volume_source, 
					   t_length, &src_tslice, 
					   &src_color_ind, &src_parity_ind, 
					   &src_corner_ind, src_seperation, 
					   j_decay);


    // Compute the solution vector for the particular source
    SystemSolverResults_t res = (*qprop)(psi, q_source);
      
    push(xml_out,"Qprop_noise");
    write(xml_out, "Noise_number" , i);
    write(xml_out, "RsdCG" , RsdCG);
    write(xml_out, "n_count", res.n_count);
    write(xml_out, "Seed" , seed);
    pop(xml_out);


    scalar_one_loop.compute(q_source,psi,i) ;
    scalar_one_kilcup_loop.compute(psi,i,Mass);
    scalar_two_loop.compute(q_source,psi,i) ;
    scalar_four_loop.compute(q_source,psi,i);
    scalar_four_kilcup_loop.compute(psi,i,Mass);
    eta3_loop.compute(q_source,psi,i) ;
    eta4_loop.compute(q_source,psi,i) ;
    eta4_kilcup_loop.compute(psi,i, Mass);
    eta0_loop.compute(q_source,psi,i) ;

     if(loop_checkpoint){
      //write each measurement to the XML file

      scalar_one_loop.dump(xml_out,i) ;
      scalar_one_kilcup_loop.dump(xml_out,i);
      scalar_two_loop.dump(xml_out,i) ;
      scalar_four_loop.dump(xml_out,i) ;
      scalar_four_kilcup_loop.dump(xml_out,i) ;
      eta3_loop.dump(xml_out,i) ;
      eta4_loop.dump(xml_out,i) ;
      eta4_kilcup_loop.dump(xml_out,i) ;
      eta0_loop.dump(xml_out,i) ;
    }

 } // Nsamples


  // write output from the loop calc
  scalar_one_loop.dump(xml_out) ;
  scalar_one_kilcup_loop.dump(xml_out);
  scalar_two_loop.dump(xml_out) ;
      scalar_four_loop.dump(xml_out) ;
      scalar_four_kilcup_loop.dump(xml_out) ;
  eta3_loop.dump(xml_out) ;
  eta4_loop.dump(xml_out) ;
  eta4_kilcup_loop.dump(xml_out) ;
  eta0_loop.dump(xml_out) ;

  // end of this section
  pop(xml_out);

}

  /**********************************************************************/


  //
  //  version used in test code.
  //

  void ks_local_loops(
		 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
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
		 VolSrc_type volume_source,
		 int src_seperation,
		 int j_decay){

  int src_tslice=0;
  int src_color_ind = 0;
  int src_parity_ind = 0;
  int src_corner_ind =0;
  Real coverage_fraction;

  push(xml_out,"ks_local_loops");

  // write common parameters
  write(xml_out, "Mass" , Mass);

  write_out_source_type(xml_out, volume_source);

  write(xml_out, "Number_of_samples" , Nsamp);


  //  parse input files


  // the wrapped disconnected loops
  bool gauge_shift ;
  bool sym_shift ;
  bool loop_checkpoint;

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
  try{
    read(xml_in, "/propagator/param/loop_checkpoint", loop_checkpoint ) ;
  }catch (const string& e){
    QDPIO::cerr << "Error reading data: " << e << endl;
    throw;
  }


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

  // set up the loop code
  local_scalar_loop                 scalar_one_loop(t_length,Nsamp,
						    u,type_of_shift) ;
  local_scalar_kilcup_loop          scalar_one_kilcup_loop(t_length, Nsamp,
							   u, type_of_shift);
  non_local_scalar_loop             scalar_two_loop(t_length,Nsamp,
						    u,type_of_shift) ;
  fourlink_scalar_loop scalar_four_loop(t_length, Nsamp,
					u, type_of_shift);
  fourlink_scalar_kilcup_loop scalar_four_kilcup_loop(t_length, Nsamp,
						      u, type_of_shift);
  threelink_pseudoscalar_loop       eta3_loop(t_length,Nsamp,
					      u,type_of_shift) ;
  fourlink_pseudoscalar_loop        eta4_loop(t_length,Nsamp,
					      u,type_of_shift) ;
  fourlink_pseudoscalar_kilcup_loop eta4_kilcup_loop(t_length, Nsamp,
						     u,type_of_shift) ;

  // for test purposes
  zerolink_pseudoscalar_loop        eta0_loop(t_length, Nsamp,
					      u,type_of_shift) ;

    // Seed the RNG with the cfg number for now
  QDP::Seed seed;
  seed = CFGNO;
  RNG::setrn(seed);

  for(int i = 0; i < Nsamp; ++i){
    psi = zero;   // note this is ``zero'' and not 0
    RNG::savern(seed);
    QDPIO::cout << "SEED = " << seed << endl;

    QDPIO::cout << "Noise sample: " << i << endl;

    // Fill the volume with random noise 
 
    coverage_fraction = fill_volume_source(q_source, volume_source, 
					   t_length, &src_tslice, 
					   &src_color_ind, &src_parity_ind, 
					   &src_corner_ind, src_seperation, 
					   j_decay);

    // Compute the solution vector for the particular source
    SystemSolverResults_t res = (*qprop)(psi, q_source);
      
    push(xml_out,"Qprop_noise");
    write(xml_out, "Noise_number" , i);
    write(xml_out, "RsdCG" , RsdCG);
    write(xml_out, "n_count", res.n_count);
    write(xml_out, "Seed" , seed);
    pop(xml_out);


    scalar_one_loop.compute(q_source,psi,i) ;
    scalar_one_kilcup_loop.compute(psi, i, Mass);
    scalar_two_loop.compute(q_source,psi,i) ;
    scalar_four_loop.compute(q_source,psi,i);
    scalar_four_kilcup_loop.compute(psi,i,Mass);
    eta3_loop.compute(q_source,psi,i) ;
    eta4_loop.compute(q_source,psi,i) ;
    eta4_kilcup_loop.compute(psi,i, Mass);
    eta0_loop.compute(q_source,psi,i) ;

    if(loop_checkpoint){
      //write each measurement to the XML file

      scalar_one_loop.dump(xml_out,i) ;
      scalar_one_kilcup_loop.dump(xml_out,i);
      scalar_two_loop.dump(xml_out,i) ;
      scalar_four_loop.dump(xml_out,i) ;
      scalar_four_kilcup_loop.dump(xml_out,i) ;
      eta3_loop.dump(xml_out,i) ;
      eta4_loop.dump(xml_out,i) ;
      eta4_kilcup_loop.dump(xml_out,i) ;
      eta0_loop.dump(xml_out,i) ;
    }

  } // Nsamples


  // write output from the loop calc
  scalar_one_loop.dump(xml_out) ;
  scalar_one_kilcup_loop.dump(xml_out);
  scalar_two_loop.dump(xml_out) ;
      scalar_four_loop.dump(xml_out) ;
      scalar_four_kilcup_loop.dump(xml_out) ;
  eta3_loop.dump(xml_out) ;
  eta4_loop.dump(xml_out) ;
  eta4_kilcup_loop.dump(xml_out) ;
  eta0_loop.dump(xml_out) ;

  // end of this section
  pop(xml_out);

}

  /**********************************************************************/

//  fuzz the loops

void ks_fuzz_loops_X(
		 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
		 LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi ,
		 LatticeStaggeredFermion & psi_fuzz ,
		 const multi1d<LatticeColorMatrix> & u,
		 const multi1d<LatticeColorMatrix> & u_smr,
		 XMLWriter & xml_out, 
		 bool gauge_shift,
		 bool sym_shift,
		 bool loop_checkpoint,
		 int t_length,
		 Real Mass,
		 int Nsamp,
		 Real RsdCG,
		 int CFGNO,
		 VolSrc_type volume_source,
		 int fuzz_width, 
		 int src_seperation,
		 int j_decay){

  int src_tslice=0;
  int src_color_ind = 0;
  int src_parity_ind = 0;
  int src_corner_ind =0;
  Real coverage_fraction;

    push(xml_out,"fuzz_loops_s");

    // write common parameters
    write(xml_out, "Mass" , Mass);


    write_out_source_type(xml_out, volume_source);


    write(xml_out, "Number_of_samples" , Nsamp);
    write(xml_out, "fuzz_width" , fuzz_width);


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



   // set up the loop code
    local_scalar_loop               scalar_one_loop(t_length,Nsamp,
						    u,type_of_shift) ;
    local_scalar_kilcup_loop        scalar_one_kilcup_loop(t_length,Nsamp,
							   u, type_of_shift);
    non_local_scalar_loop           scalar_two_loop(t_length,Nsamp,
						    u,type_of_shift) ;
    fourlink_scalar_loop scalar_four_loop(t_length, Nsamp,
					  u, type_of_shift);
    fourlink_scalar_kilcup_loop scalar_four_kilcup_loop(t_length, Nsamp,
							u, type_of_shift);
    threelink_pseudoscalar_loop     eta3_loop(t_length,Nsamp,
					      u,type_of_shift) ;
    fourlink_pseudoscalar_loop      eta4_loop(t_length,Nsamp,
					      u,type_of_shift) ;

  fourlink_pseudoscalar_kilcup_loop eta4_kilcup_loop(t_length, Nsamp,
						     u,type_of_shift) ;

  // for test purposes
  zerolink_pseudoscalar_loop        eta0_loop(t_length, Nsamp,
					      u,type_of_shift) ;

    // Seed the RNG with the cfg number for now
    QDP::Seed seed;
    seed = CFGNO;
    RNG::setrn(seed);


  for(int i = 0; i < Nsamp; ++i){
    psi = zero;   // note this is ``zero'' and not 0
    RNG::savern(seed);
    QDPIO::cout << "SEED = " << seed << endl;

    QDPIO::cout << "Noise sample: " << i << endl;

    // Fill the volume with random noise 
    coverage_fraction = fill_volume_source(q_source, volume_source, 
					   t_length, &src_tslice, 
					   &src_color_ind, &src_parity_ind, 
					   &src_corner_ind, src_seperation, 
					   j_decay);

    // Compute the solution vector for the particular source
    SystemSolverResults_t res = (*qprop)(psi, q_source);
      
    push(xml_out,"Qprop_noise");
    write(xml_out, "Noise_number" , i);
    write(xml_out, "RsdCG" , RsdCG);
    write(xml_out, "n_count", res.n_count);
    write(xml_out, "Seed" , seed);
    pop(xml_out);


    fuzz_smear(u_smr, psi,psi_fuzz,
	       fuzz_width, j_decay) ;


    scalar_one_loop.compute(q_source,psi_fuzz,i) ;
    scalar_one_kilcup_loop.compute(psi, i, Mass);
    scalar_two_loop.compute(q_source,psi_fuzz,i) ;
    scalar_four_loop.compute(q_source,psi,i);
    scalar_four_kilcup_loop.compute(psi,i,Mass);
    eta3_loop.compute(q_source,psi_fuzz,i) ;
    eta4_loop.compute(q_source,psi_fuzz,i) ;
    eta4_kilcup_loop.compute(psi,i, Mass);
    eta0_loop.compute(q_source,psi,i) ;


    if(loop_checkpoint){
      //write each measurement to the XML file

      scalar_one_loop.dump(xml_out,i) ;
      scalar_one_kilcup_loop.dump(xml_out,i);
      scalar_two_loop.dump(xml_out,i) ;
      scalar_four_loop.dump(xml_out,i) ;
      scalar_four_kilcup_loop.dump(xml_out,i) ;
      eta3_loop.dump(xml_out,i) ;
      eta4_loop.dump(xml_out,i) ;
      eta4_kilcup_loop.dump(xml_out,i) ;
      eta0_loop.dump(xml_out,i) ;
    }

  } // Nsamples


  // write output from the loop calc
  scalar_one_loop.dump(xml_out) ;
  scalar_one_kilcup_loop.dump(xml_out);
  scalar_two_loop.dump(xml_out) ;
      scalar_four_loop.dump(xml_out) ;
      scalar_four_kilcup_loop.dump(xml_out) ;
  eta3_loop.dump(xml_out) ;
  eta4_loop.dump(xml_out) ;
  eta4_kilcup_loop.dump(xml_out) ;
  eta0_loop.dump(xml_out) ;

  // end of this section
  pop(xml_out);

}

  /**********************************************************************/

//  fuzz the loops
//  HACK

void ks_fuzz_loops(
		 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
		 LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi ,
		 LatticeStaggeredFermion & psi_fuzz ,
		 const multi1d<LatticeColorMatrix> & u,
		 const multi1d<LatticeColorMatrix> & u_smr,
		 XMLWriter & xml_out, 
		 bool gauge_shift,
		 bool sym_shift,
		 bool loop_checkpoint,
		 int t_length,
		 Real Mass,
		 int Nsamp,
		 Real RsdCG,
		 int CFGNO,
		 VolSrc_type volume_source,
		 int fuzz_width, 
		 int src_seperation,
		 int j_decay, bool binary_loop_checkpoint,
                 std::string binary_name){

  int src_tslice=0;
  int src_color_ind = 0;
  int src_parity_ind = 0;
  int src_corner_ind =0;
  Real coverage_fraction;
  int j;
  int fuzz_index;

    push(xml_out,"fuzz_loops_s");

    // write common parameters
    write(xml_out, "Mass" , Mass);


    write_out_source_type(xml_out, volume_source);


    write(xml_out, "Number_of_samples" , Nsamp);
    write(xml_out, "fuzz_width" , fuzz_width);


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

   // set up the loop code
    local_scalar_loop               scalar_one_loop(t_length,Nsamp,
						    u,type_of_shift) ;
    local_scalar_kilcup_loop        scalar_one_kilcup_loop(t_length, Nsamp,
							   u, type_of_shift);
    non_local_scalar_loop           scalar_two_loop(t_length,Nsamp,
						    u,type_of_shift) ;
    fourlink_scalar_loop scalar_four_loop(t_length, Nsamp,
					  u, type_of_shift);
    fourlink_scalar_kilcup_loop scalar_four_kilcup_loop(t_length, Nsamp,
							u, type_of_shift);
    threelink_pseudoscalar_loop     eta3_loop(t_length,Nsamp,
					      u,type_of_shift) ;
    fourlink_pseudoscalar_loop      eta4_loop(t_length,Nsamp,
					      u,type_of_shift) ;

    fourlink_pseudoscalar_kilcup_loop eta4_kilcup_loop(t_length, Nsamp,
						       u,type_of_shift) ;

    // for test purposes
    zerolink_pseudoscalar_loop        eta0_loop(t_length, Nsamp,
						u,type_of_shift);



    // fuzzed loops here
    local_scalar_loop_fuzz       scalar_one_loop_fuzz(t_length, Nsamp, u,
						      type_of_shift);
    local_scalar_kilcup_loop_fuzz scalar_one_kilcup_loop_fuzz(t_length, Nsamp,
							      u, type_of_shift);
    non_local_scalar_loop_fuzz   scalar_two_loop_fuzz(t_length,Nsamp, u,
						      type_of_shift);

    threelink_pseudoscalar_loop_fuzz     eta3_loop_fuzz(t_length,Nsamp,
							u,type_of_shift);
    fourlink_pseudoscalar_loop_fuzz      eta4_loop_fuzz(t_length,Nsamp,
							u,type_of_shift);

    fourlink_pseudoscalar_kilcup_loop_fuzz eta4_kilcup_loop_fuzz(t_length,
								 Nsamp, u,
								type_of_shift);



    // set up the loop code

    // Seed the RNG with the cfg number for now
    QDP::Seed seed;
    seed = CFGNO;
    RNG::setrn(seed);


    for(int i = 0; i < Nsamp; ++i){
      psi = zero;                // note this is ``zero'' and not 0
      RNG::savern(seed);
      QDPIO::cout << "SEED = " << seed << endl;

      QDPIO::cout << "Noise sample: " << i << endl;

      // Fill the volume with random noise 
      coverage_fraction = fill_volume_source(q_source, volume_source,
					     t_length, &src_tslice,
					     &src_color_ind, &src_parity_ind,
					     &src_corner_ind, src_seperation,
					     j_decay);
      SystemSolverResults_t res  ;
      {
        StopWatch swatch;
        swatch.start();

      // Compute the solution vector for the particular source
      res = (*qprop)(psi, q_source);
      
        swatch.stop();
        double time_in_sec  = swatch.getTimeInSeconds();
	QDPIO::cout << "ks_fuzz_loops::PROF INVERTER  [" << i << "] " << time_in_sec << " sec" << endl;



      }


      push(xml_out,"Qprop_noise");
      write(xml_out, "Noise_number" , i);
      write(xml_out, "RsdCG" , RsdCG);
      write(xml_out, "n_count", res.n_count);
      write(xml_out, "Seed" , seed);
      pop(xml_out);

      {
        StopWatch swatch;
        swatch.start();

      fuzz_smear(u_smr, psi,psi_fuzz, fuzz_width, j_decay) ;

        swatch.stop();
        double time_in_sec  = swatch.getTimeInSeconds();
	QDPIO::cout << "ks_fuzz_loops::PROF fuzz_smear  [" << i << "] " << time_in_sec << " sec" << endl;


      }

      {
        StopWatch swatch;
        swatch.start();


      // compute the un-fuzzed operators
      scalar_one_loop.compute(q_source,psi,i);
      scalar_one_kilcup_loop.compute(psi, i, Mass);
      scalar_two_loop.compute(q_source,psi,i);
      scalar_four_loop.compute(q_source,psi,i);
      scalar_four_kilcup_loop.compute(psi,i,Mass);
      eta3_loop.compute(q_source,psi,i);
      eta4_loop.compute(q_source,psi,i);
      eta4_kilcup_loop.compute(psi,i, Mass);
      eta0_loop.compute(q_source,psi,i);

      // now compute the fuzzed operators
      scalar_one_loop_fuzz.compute(q_source,psi_fuzz,i);
      scalar_one_kilcup_loop_fuzz.compute(psi_fuzz, psi, i, Mass);
      scalar_two_loop_fuzz.compute(q_source,psi_fuzz,i) ;
      eta3_loop_fuzz.compute(q_source,psi_fuzz,i) ;
      eta4_loop_fuzz.compute(q_source,psi_fuzz,i) ;
      eta4_kilcup_loop_fuzz.compute(psi_fuzz,psi,i, Mass);


        swatch.stop();
        double time_in_sec  = swatch.getTimeInSeconds();
	QDPIO::cout << "ks_fuzz_loops::PROF compute  [" << i << "] " << time_in_sec << " sec" << endl;


      }


      if(loop_checkpoint){
        StopWatch swatch;
        swatch.start();

        //write each measurement to the XML file

	scalar_one_loop.dump(xml_out,i) ;
	scalar_one_kilcup_loop.dump(xml_out,i);
	scalar_two_loop.dump(xml_out,i) ;
      scalar_four_loop.dump(xml_out,i) ;
      scalar_four_kilcup_loop.dump(xml_out,i) ;
	eta3_loop.dump(xml_out,i) ;
	eta4_loop.dump(xml_out,i) ;
	eta4_kilcup_loop.dump(xml_out,i) ;
        eta0_loop.dump(xml_out,i) ;

	scalar_one_loop_fuzz.dump(xml_out,i) ;
	scalar_one_kilcup_loop_fuzz.dump(xml_out,i);
	scalar_two_loop_fuzz.dump(xml_out,i) ;
	eta3_loop_fuzz.dump(xml_out,i) ;
	eta4_loop_fuzz.dump(xml_out,i) ;
	eta4_kilcup_loop_fuzz.dump(xml_out,i) ;

        swatch.stop();
        double time_in_sec  = swatch.getTimeInSeconds();
	QDPIO::cout << "ks_fuzz_loops::PROF CHECKPOINT  [" << i << "] " << time_in_sec << " sec" << endl;


      }
    
    
    } // Nsamples

    // write output from the loop calc
    {
      StopWatch swatch;
      swatch.start();

  scalar_one_loop.dump(xml_out) ;
  scalar_one_kilcup_loop.dump(xml_out);
  scalar_two_loop.dump(xml_out) ;
  scalar_four_loop.dump(xml_out) ;
  scalar_four_kilcup_loop.dump(xml_out) ;
  eta3_loop.dump(xml_out) ;
  eta4_loop.dump(xml_out) ;
  eta4_kilcup_loop.dump(xml_out) ;

  eta0_loop.dump(xml_out) ;

  scalar_one_loop_fuzz.dump(xml_out);
  scalar_one_kilcup_loop_fuzz.dump(xml_out);
  scalar_two_loop_fuzz.dump(xml_out) ;
  eta3_loop_fuzz.dump(xml_out) ;
  eta4_loop_fuzz.dump(xml_out) ;
  eta4_kilcup_loop_fuzz.dump(xml_out) ;


  swatch.stop();
  double time_in_sec  = swatch.getTimeInSeconds();
  QDPIO::cout << "ks_fuzz_loops::FINAL IO  " << time_in_sec << " sec" << endl;

    }

    if(binary_loop_checkpoint )
    {
      // BINARY DUMP 
      StopWatch swatch;
      swatch.start();

      scalar_one_loop.binary_dump(binary_name) ;
      scalar_one_kilcup_loop.binary_dump(binary_name) ;
      scalar_two_loop.binary_dump(binary_name) ;
      scalar_four_loop.binary_dump(binary_name);
      scalar_four_kilcup_loop.binary_dump(binary_name);
      eta3_loop.binary_dump(binary_name) ;
      eta4_loop.binary_dump(binary_name) ;
      eta4_kilcup_loop.binary_dump(binary_name) ;
      eta0_loop.binary_dump(binary_name) ;

      scalar_one_loop_fuzz.binary_dump(binary_name) ;
      scalar_one_kilcup_loop_fuzz.binary_dump(binary_name) ;
      scalar_two_loop_fuzz.binary_dump(binary_name) ;
      eta3_loop_fuzz.binary_dump(binary_name) ;
      eta4_loop_fuzz.binary_dump(binary_name) ;
      eta4_kilcup_loop_fuzz.binary_dump(binary_name) ;

      swatch.stop();
      double time_in_sec  = swatch.getTimeInSeconds();
      QDPIO::cout << "ks_fuzz_loops::BINARY IO  " << time_in_sec << " sec" << endl;

    }


  // end of this section
  pop(xml_out);

}


/**********************************************************************/


void ks_local_loops_and_stoch_conn(
		 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
		 LatticeStaggeredFermion & q_source1, 
		 LatticeStaggeredFermion & psi1 ,
		 const multi1d<LatticeColorMatrix> & u,
		 XMLWriter & xml_out, 
		 bool gauge_shift,
		 bool sym_shift,
		 bool loop_checkpoint,
		 int t_length,
		 Real Mass,
		 int Nsamp,
		 Real RsdCG,
		 int CFGNO,
		 VolSrc_type volume_source,
		 int src_seperation,
		 int j_decay){


  LatticeStaggeredFermion psi2 ;
  LatticeStaggeredFermion q_source2 ;


  push(xml_out,"local_loops_s");

  // write common parameters
  write(xml_out, "Mass" , Mass);

  write_out_source_type(xml_out, volume_source);

  write(xml_out, "Number_of_samples" , Nsamp);

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

  // set up the loop code
  local_scalar_loop                 scalar_one_loop(t_length,Nsamp,
						    u,type_of_shift) ;
  local_scalar_kilcup_loop          scalar_one_kilcup_loop(t_length, Nsamp,
							   u, type_of_shift);
  non_local_scalar_loop             scalar_two_loop(t_length,Nsamp,
						    u,type_of_shift) ;
    fourlink_scalar_loop scalar_four_loop(t_length, Nsamp,
					  u, type_of_shift);
    fourlink_scalar_kilcup_loop scalar_four_kilcup_loop(t_length, Nsamp,
							u, type_of_shift);
  threelink_pseudoscalar_loop       eta3_loop(t_length,Nsamp,
					      u,type_of_shift) ;
  fourlink_pseudoscalar_loop        eta4_loop(t_length,Nsamp,
					      u,type_of_shift) ;

  fourlink_pseudoscalar_kilcup_loop eta4_kilcup_loop(t_length,Nsamp,
						     u,type_of_shift) ;

  zerolink_pseudoscalar_loop        eta0_loop(t_length, Nsamp,
					      u,type_of_shift) ;

  fourlink_pseudoscalar_stoch_conn   eta4_conn(t_length,Nsamp,
					       u,type_of_shift) ;

  // Seed the RNG with the cfg number for now
  QDP::Seed seed;
  seed = CFGNO;
  RNG::setrn(seed);

  int src_tslice=0;
  int src_color_ind = 0;
  int src_parity_ind = 0;
  int src_corner_ind =0;

  int src_tslice2=0;
  int src_color_ind2 = 0;
  int src_parity_ind2 = 0;
  int src_corner_ind2 =0;

  Real coverage_fraction;

  for(int i = 0; i < Nsamp; ++i){
    psi1 = zero;   // note this is ``zero'' and not 0
    psi2 = zero;   // note this is ``zero'' and not 0
    RNG::savern(seed);
    QDPIO::cout << "SEED = " << seed << endl;

    QDPIO::cout << "Noise sample: " << i << endl;

    // Fill the volume with random noise 
    coverage_fraction = fill_volume_source(q_source1, volume_source, 
					   t_length, &src_tslice, 
					   &src_color_ind, &src_parity_ind, 
					   &src_corner_ind, src_seperation, 
					   j_decay);

    //fill 2nd source for stochastic connected correlators
    coverage_fraction = fill_volume_source(q_source2, volume_source, 
					   t_length, &src_tslice2, 
					   &src_color_ind2, &src_parity_ind2, 
					   &src_corner_ind2, src_seperation, 
					   j_decay);

    // Compute the solution vector for the particular source
    SystemSolverResults_t res1 = (*qprop)(psi1, q_source1);
    SystemSolverResults_t res2 = (*qprop)(psi2, q_source2);
      
    push(xml_out,"Qprop_noise");
    write(xml_out, "Noise_number" , i);
    write(xml_out, "RsdCG" , RsdCG);
    write(xml_out, "n_count1", res1.n_count);
    write(xml_out, "n_count2", res2.n_count);
    write(xml_out, "Seed" , seed);
    pop(xml_out);


    scalar_one_loop.compute(q_source1,psi1,i) ;
    scalar_one_kilcup_loop.compute(psi1,i, Mass);  
    scalar_two_loop.compute(q_source1,psi1,i) ;
    scalar_four_loop.compute(q_source1,psi1,i);
    scalar_four_kilcup_loop.compute(psi1,i,Mass);
    eta3_loop.compute(q_source1,psi1,i) ;
    eta4_loop.compute(q_source1,psi1,i) ;
    eta4_kilcup_loop.compute(psi1,i, Mass);
    eta0_loop.compute(q_source1,psi1,i) ;

    eta4_conn.compute(q_source1, q_source2, psi1, psi2, i) ;

    if(loop_checkpoint){
      //write each measurement to the XML file

      scalar_one_loop.dump(xml_out,i) ;
      scalar_one_kilcup_loop.dump(xml_out,i);
      scalar_two_loop.dump(xml_out,i) ;
      scalar_four_loop.dump(xml_out,i) ;
      scalar_four_kilcup_loop.dump(xml_out,i) ;
      eta3_loop.dump(xml_out,i) ;
      eta4_loop.dump(xml_out,i) ;
      eta4_kilcup_loop.dump(xml_out,i) ;
      eta0_loop.dump(xml_out,i) ;
    }

    // MUST checkpoint stochastic connected correlator measurements!
    eta4_conn.dump(xml_out,i) ;
    printf("OUTNOW!!!!\n");fflush(stdout);
 } // Nsamples


  // write output from the loop calc
  scalar_one_loop.dump(xml_out) ;
  scalar_one_kilcup_loop.dump(xml_out);
  scalar_two_loop.dump(xml_out) ;
      scalar_four_loop.dump(xml_out) ;
      scalar_four_kilcup_loop.dump(xml_out) ;
  eta3_loop.dump(xml_out) ;
  eta4_loop.dump(xml_out) ;
  eta4_kilcup_loop.dump(xml_out) ;
  eta0_loop.dump(xml_out) ;

  // no point in dumping connected correlator measurements, 
  // since corrs must be formed with each noise source, but what the hell, 
  // maybe there is a diagnostic use for it.
    printf("OUTNOW1!!!!\n");fflush(stdout);
  eta4_conn.dump(xml_out) ;

    printf("OUTNOW2!!!!\n");fflush(stdout);
  // end of this section
  pop(xml_out);

}

  /**********************************************************************/

  /**********************************************************************/

//  fuzz the loops

void ks_fuzz_loops_stoch_conn(
		 Handle< SystemSolver<LatticeStaggeredFermion> > & qprop,
		 LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi ,
		 LatticeStaggeredFermion & psi_fuzz ,
		 const multi1d<LatticeColorMatrix> & u,
		 const multi1d<LatticeColorMatrix> & u_smr,
		 XMLWriter & xml_out, 
		 bool gauge_shift,
		 bool sym_shift,
		 bool loop_checkpoint,
		 int t_length,
		 Real Mass,
		 int Nsamp,
		 Real RsdCG,
		 int CFGNO,
		 VolSrc_type volume_source,
		 int fuzz_width, 
		 int src_seperation,
		 int j_decay){

  int src_tslice=0;
  int src_color_ind = 0;
  int src_parity_ind = 0;
  int src_corner_ind =0;

  int src_tslice2=0;
  int src_color_ind2 = 0;
  int src_parity_ind2 = 0;
  int src_corner_ind2 =0;

  Real coverage_fraction;
  int j;
  int fuzz_index;

  LatticeStaggeredFermion psi2 ;
  LatticeStaggeredFermion q_source2 ;

  bool fuzz_sink=false;
  bool fuzz_src=false;

    push(xml_out,"fuzz_loops_s");

    // write common parameters
    write(xml_out, "Mass" , Mass);


    write_out_source_type(xml_out, volume_source);


    write(xml_out, "Number_of_samples" , Nsamp);
    write(xml_out, "fuzz_width" , fuzz_width);


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

   // set up the loop code
    local_scalar_loop               scalar_one_loop(t_length,Nsamp,
						    u,type_of_shift) ;
    local_scalar_kilcup_loop        scalar_one_kilcup_loop(t_length, Nsamp,
							   u,type_of_shift);
    non_local_scalar_loop           scalar_two_loop(t_length,Nsamp,
						    u,type_of_shift) ;
    fourlink_scalar_loop scalar_four_loop(t_length, Nsamp,
					  u, type_of_shift);
    fourlink_scalar_kilcup_loop scalar_four_kilcup_loop(t_length, Nsamp,
							u, type_of_shift);
    threelink_pseudoscalar_loop     eta3_loop(t_length,Nsamp,
					      u,type_of_shift) ;
    fourlink_pseudoscalar_loop      eta4_loop(t_length,Nsamp,
					      u,type_of_shift) ;

    fourlink_pseudoscalar_kilcup_loop eta4_kilcup_loop(t_length, Nsamp,
						       u,type_of_shift) ;

    // for test purposes
    zerolink_pseudoscalar_loop        eta0_loop(t_length, Nsamp,
						u,type_of_shift);



    // fuzzed loops here
    local_scalar_loop_fuzz       scalar_one_loop_fuzz(t_length, Nsamp, u,
						      type_of_shift);
    local_scalar_kilcup_loop_fuzz scalar_one_kilcup_loop_fuzz(t_length, Nsamp,
							      u, type_of_shift);
    non_local_scalar_loop_fuzz   scalar_two_loop_fuzz(t_length,Nsamp, u,
						      type_of_shift);
    threelink_pseudoscalar_loop_fuzz     eta3_loop_fuzz(t_length,Nsamp,
							u,type_of_shift);
    fourlink_pseudoscalar_loop_fuzz      eta4_loop_fuzz(t_length,Nsamp,
							u,type_of_shift);

    fourlink_pseudoscalar_kilcup_loop_fuzz eta4_kilcup_loop_fuzz(t_length,
								 Nsamp, u,
								type_of_shift);


    //  fourlink_pseudoscalar_stoch_conn   eta4_conn(t_length,Nsamp,
    //					       u,type_of_shift) ;


  fourlink_pseudoscalar_stoch_conn   eta4_conn(t_length,Nsamp,
					       u,type_of_shift,
					       false,false) ;

  fourlink_pseudoscalar_stoch_conn   eta4_conn_FsrcLsink(t_length,Nsamp,
							 u,type_of_shift,
							 true,false) ;

  fourlink_pseudoscalar_stoch_conn   eta4_conn_LsrcFsink(t_length,Nsamp,
							 u,type_of_shift, 
							 false,true) ;

    // set up the loop code

    // Seed the RNG with the cfg number for now
    QDP::Seed seed;
    seed = CFGNO;
    RNG::setrn(seed);


    for(int i = 0; i < Nsamp; ++i){
      psi = zero;                // note this is ``zero'' and not 0
      psi2 = zero;   // note this is ``zero'' and not 0

      RNG::savern(seed);
      QDPIO::cout << "SEED = " << seed << endl;

      QDPIO::cout << "Noise sample: " << i << endl;

      // Fill the volume with random noise 
      coverage_fraction = fill_volume_source(q_source, volume_source,
					     t_length, &src_tslice,
					     &src_color_ind, &src_parity_ind,
					     &src_corner_ind, src_seperation,
					     j_decay);

      //fill 2nd source for stochastic connected correlators
      coverage_fraction = fill_volume_source(q_source2, volume_source, 
					     t_length, &src_tslice2,
					     &src_color_ind2, &src_parity_ind2,
					     &src_corner_ind2, src_seperation,
					     j_decay);


      // Compute the solution vector for the particular source
      SystemSolverResults_t res = (*qprop)(psi, q_source);
      SystemSolverResults_t res2 = (*qprop)(psi2, q_source2);

      push(xml_out,"Qprop_noise");
      write(xml_out, "Noise_number" , i);
      write(xml_out, "RsdCG" , RsdCG);
      write(xml_out, "n_count", res.n_count);
      write(xml_out, "Seed" , seed);
      pop(xml_out);

      fuzz_smear(u_smr, psi,psi_fuzz, fuzz_width, j_decay) ;


      // compute the un-fuzzed operators
      scalar_one_loop.compute(q_source,psi,i);
      scalar_one_kilcup_loop.compute(psi, i, Mass);
      scalar_two_loop.compute(q_source,psi,i);
      scalar_four_loop.compute(q_source,psi,i);
      scalar_four_kilcup_loop.compute(psi,i,Mass);
      eta3_loop.compute(q_source,psi,i);
      eta4_loop.compute(q_source,psi,i);
      eta4_kilcup_loop.compute(psi,i, Mass);
      eta0_loop.compute(q_source,psi,i);

      //      eta4_conn.compute(q_source, q_source2, psi, psi2, i) ;

      eta4_conn.compute(q_source, q_source2,
			psi, psi2,
			u_smr, fuzz_width, j_decay, i);

      eta4_conn_FsrcLsink.compute(q_source, q_source2,
				  psi, psi2,
				  u_smr, fuzz_width, j_decay, i);
      eta4_conn_LsrcFsink.compute(q_source, q_source2,
				  psi, psi2,
				  u_smr, fuzz_width, j_decay, i);

      // now compute the fuzzed operators
      scalar_one_loop_fuzz.compute(q_source,psi_fuzz,i);
      scalar_one_kilcup_loop_fuzz.compute(psi_fuzz, psi, i, Mass);
      scalar_two_loop_fuzz.compute(q_source,psi_fuzz,i) ;
      eta3_loop_fuzz.compute(q_source,psi_fuzz,i) ;
      eta4_loop_fuzz.compute(q_source,psi_fuzz,i) ;
      eta4_kilcup_loop_fuzz.compute(psi_fuzz,psi,i, Mass);




      if(loop_checkpoint){
        //write each measurement to the XML file

	scalar_one_loop.dump(xml_out,i) ;
	scalar_one_kilcup_loop.dump(xml_out,i);
	scalar_two_loop.dump(xml_out,i) ;
	scalar_four_loop.dump(xml_out,i) ;
	scalar_four_kilcup_loop.dump(xml_out,i); 
	eta3_loop.dump(xml_out,i) ;
	eta4_loop.dump(xml_out,i) ;
	eta4_kilcup_loop.dump(xml_out,i) ;
        eta0_loop.dump(xml_out,i) ;

	scalar_one_loop_fuzz.dump(xml_out,i) ;
	scalar_one_kilcup_loop_fuzz.dump(xml_out,i);
	scalar_two_loop_fuzz.dump(xml_out,i) ;
	eta3_loop_fuzz.dump(xml_out,i) ;
	eta4_loop_fuzz.dump(xml_out,i) ;
	eta4_kilcup_loop_fuzz.dump(xml_out,i) ;
      }
    
      // MUST checkpoint stochastic connected correlator measurements!
      eta4_conn.dump(xml_out,i) ;
      eta4_conn_FsrcLsink.dump(xml_out,i) ;
      eta4_conn_LsrcFsink.dump(xml_out,i) ;

      //      printf("OUTNOW!!!!\n");fflush(stdout);
     
    } // Nsamples

    // write output from the loop calc
  scalar_one_loop.dump(xml_out) ;
  scalar_one_kilcup_loop.dump(xml_out);
  scalar_two_loop.dump(xml_out) ;
      scalar_four_loop.dump(xml_out) ;
      scalar_four_kilcup_loop.dump(xml_out) ;
  eta3_loop.dump(xml_out) ;
  eta4_loop.dump(xml_out) ;
  eta4_kilcup_loop.dump(xml_out) ;

  eta0_loop.dump(xml_out) ;

  scalar_one_loop_fuzz.dump(xml_out) ;
  scalar_one_kilcup_loop_fuzz.dump(xml_out);
  scalar_two_loop_fuzz.dump(xml_out) ;
  eta3_loop_fuzz.dump(xml_out) ;
  eta4_loop_fuzz.dump(xml_out) ;
  eta4_kilcup_loop_fuzz.dump(xml_out) ;

  // no point in dumping connected correlator measurements, 
  // since corrs must be formed with each noise source, but what the hell, 
  // maybe there is a diagnostic use for it.
    printf("OUTNOW1!!!!\n");fflush(stdout);
    eta4_conn.dump(xml_out) ;
    eta4_conn_FsrcLsink.dump(xml_out) ;
    eta4_conn_LsrcFsink.dump(xml_out) ;

  // end of this section
  pop(xml_out);

}


/**********************************************************************/


}  // end namespace Chroma
