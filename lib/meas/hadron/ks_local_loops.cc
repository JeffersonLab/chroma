/* + */
/* $Id: ks_local_loops.cc,v 2.1 2005-11-01 04:15:16 edwards Exp $ ($Date: 2005-11-01 04:15:16 $) */


#include "fermact.h"
#include "meas/hadron/ks_local_loops.h"
#include "meas/hadron/hadron_s.h"
#include "meas/sources/z2_src.h"

#include "meas/smear/fuzz_smear.h"

namespace Chroma {



void ks_local_loops(
		 Handle<const SystemSolver<LatticeStaggeredFermion> > & qprop,
		 LatticeStaggeredFermion & q_source, 
		 LatticeStaggeredFermion & psi ,
		 const multi1d<LatticeColorMatrix> & u,
		 XMLWriter & xml_out, 
		 bool gauge_shift,
		 int t_length,
		 Real Mass,
		 int Nsamp,
		 Real RsdCG,
		 int CFGNO,
		 int volume_source
		 )
{

    push(xml_out,"local_loops_s");

    // write common parameters
    write(xml_out, "Mass" , Mass);
    if( volume_source == Z2NOISE  )
      write(xml_out, "Random_volume_source" , "Z2NOISE");
    else if( volume_source == GAUSSIAN )
      write(xml_out, "Random_volume_source" , "GAUSSIAN");

    write(xml_out, "Number_of_samples" , Nsamp);


    Stag_shift_option type_of_shift ; 
    if( gauge_shift )
      type_of_shift = SYM_GAUGE_INVAR ; 
    else 
      type_of_shift = NON_GAUGE_INVAR ; 

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
    try
      {
	read(xml_in, "/propagator/param/use_gauge_invar_oper", gauge_shift ) ;
      }
    catch (const string& e)
      {
	QDPIO::cerr << "Error reading data: " << e << endl;
	throw;
      }


    Stag_shift_option type_of_shift ; 
    if( gauge_shift )
      type_of_shift = SYM_GAUGE_INVAR ; 
    else 
      type_of_shift = NON_GAUGE_INVAR ; 

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

    Stag_shift_option type_of_shift ; 
    if( gauge_shift )
      type_of_shift = SYM_GAUGE_INVAR ; 
    else 
      type_of_shift = NON_GAUGE_INVAR ; 

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
