// $Id: util_compute_meson_s.cc,v 3.1 2009-04-26 21:02:51 mcneile Exp $
/*! \file
 * \brief Wrapper code to compute staggered meson correlators.
 *
 * Spectrum calculations
 */

#include "handle.h"
//#include "actions/ferm/fermbcs/fermbcs.h"
#include "actions/ferm/fermacts/fermacts_s.h"
#include "meas/hadron/hadron_s.h"
#include "meas/smear/fuzz_smear.h"
#include "meas/sources/srcfil.h"
#include "util/ferm/transf.h"
#include "meas/hadron/pion_local_s.h"
#include "meas/hadron/scalar_local_s.h"

#include "util_compute_quark_prop_s.h"

namespace Chroma 
{ 


  void compute_vary_ps(
		      LatticeStaggeredPropagator & quark_propagator_Lsink_Lsrc,
		      LatticeStaggeredPropagator & quark_propagator_Fsink_Lsrc,
		      LatticeStaggeredPropagator & quark_propagator_Lsink_Fsrc,
		      LatticeStaggeredPropagator & quark_propagator_Fsink_Fsrc,
		      const multi1d<LatticeColorMatrix> & u , 
		      bool gauge_shift, bool sym_shift ,
		      XMLWriter & xml_out,
		      int j_decay,int t_length, int t_source){

    Stag_shift_option type_of_shift;

    if((gauge_shift) && (sym_shift)){
      type_of_shift=SYM_GAUGE_INVAR;
    }else{
      if((gauge_shift) && (!sym_shift)){
	type_of_shift=GAUGE_INVAR;
      }else{
	if((!gauge_shift) && (sym_shift)){
	  type_of_shift=SYM_NON_GAUGE_INVAR;
	}else{
	  if((!gauge_shift) && (!sym_shift)){
	    type_of_shift=NON_GAUGE_INVAR ;
	  }
	}
      }
    }


  staggered_local_pion pion(t_length,u, type_of_shift) ;


push(xml_out,"Meson_correlators");

  // ---------- LL ----------
  pion.compute(quark_propagator_Lsink_Lsrc,
		  quark_propagator_Lsink_Lsrc,j_decay);


  push(xml_out, "Lsink_Lsrc");

  pion.dump(t_source,xml_out);

  pop(xml_out);
  
  // ---------- FL ----------
  pion.compute(quark_propagator_Lsink_Lsrc,
		  quark_propagator_Fsink_Lsrc,j_decay);
  push(xml_out, "Fsink_Lsrc");
  pion.dump(t_source,xml_out);
  pop(xml_out);

  
  // ---------- LF ----------
  pion.compute(quark_propagator_Lsink_Lsrc,
		  quark_propagator_Lsink_Fsrc,j_decay);
  push(xml_out, "Lsink_Fsrc");
  pion.dump(t_source,xml_out);
  pop(xml_out);
  

  // ---------- FF ----------
  pion.compute(quark_propagator_Lsink_Lsrc,
		  quark_propagator_Fsink_Fsrc,j_decay);
  push(xml_out, "Fsink_Fsrc");
  pion.dump(t_source,xml_out);
  pop(xml_out);

  // ---------------------------------------
  pop(xml_out);

  QDPIO::cout << "Computed variational fuzzed mesons"  << endl;

  }

  void compute_vary_scalar(LatticeStaggeredPropagator& prop_Lsn_Lsr,
			   LatticeStaggeredPropagator& prop_Fsn_Lsr,
			   LatticeStaggeredPropagator& prop_Lsn_Fsr,
			   LatticeStaggeredPropagator& prop_Fsn_Fsr,
			   const multi1d<LatticeColorMatrix>& u,
			   bool gauge_shift, bool sym_shift,
			   XMLWriter& xml_out, int j_decay, int t_length, int t_source) 
  {
    Stag_shift_option type_of_shift;

    if (gauge_shift && sym_shift) {
      type_of_shift = SYM_GAUGE_INVAR;
    } else {
      if ( gauge_shift && !sym_shift ) {
	type_of_shift = GAUGE_INVAR;
      } else {
	if ( !gauge_shift && sym_shift ) {
	  type_of_shift = SYM_NON_GAUGE_INVAR;
	} else {
	  if ( !gauge_shift && !sym_shift ) {
	    type_of_shift = NON_GAUGE_INVAR;
	  }
	}
      }
    }
  
    staggered_local_scalar scalar(t_length, u, type_of_shift);

    push(xml_out, "Scalar_multichannel");
    scalar.compute(prop_Lsn_Lsr, prop_Lsn_Lsr, j_decay);
    push(xml_out, "Lsink_Lsrc");
    scalar.dump(t_source,xml_out);
    pop(xml_out);

    scalar.compute(prop_Lsn_Lsr, prop_Lsn_Fsr, j_decay);
    push(xml_out, "Fsink_Lsrc");
    scalar.dump(t_source,xml_out);
    pop(xml_out);

    scalar.compute(prop_Lsn_Lsr, prop_Lsn_Fsr, j_decay);
    push(xml_out, "Lsink_Fsrc");
    scalar.dump(t_source,xml_out);
    pop(xml_out);

    scalar.compute(prop_Lsn_Lsr, prop_Fsn_Fsr, j_decay);
    push(xml_out, "Fsink_Fsrc");
    scalar.dump(t_source,xml_out);
    pop(xml_out);

    pop(xml_out);

    QDPIO::cout << "Computed multichannel local scalar\n";
  }


} // end of namespace
