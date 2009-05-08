// $Id: util_compute_8_pions_s.cc,v 3.1 2009-05-08 10:15:40 mcneile Exp $
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
#include "meas/hadron/pions_s.h"
#include "meas/hadron/stag_propShift_s.h"


#include "util_compute_quark_prop_s.h"

namespace Chroma { 


  void compute_8_pions(	multi1d<LatticeStaggeredPropagator> & stag_prop,
		       const multi1d<LatticeColorMatrix> & u , 
		       bool gauge_shift, bool sym_shift,
		       XMLWriter & xml_out,
		       int j_decay, int t_length, int t_source,
		       bool binary_meson_dump, std::string binary_name){

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

    staggered_pions pseudoscalar(t_length,  u, type_of_shift) ;



  // ---------- LL ----------
    pseudoscalar.compute(stag_prop,j_decay);


    if(binary_meson_dump){

      std::string tagged_filename_base;
      tagged_filename_base=binary_name+"PS.LL.";
      pseudoscalar.binary_dump(t_source,tagged_filename_base);

    }else{

      push(xml_out,"Pseudoscalars");
      push(xml_out, "Lsink_Lsrc");
      pseudoscalar.dump(t_source,xml_out);
      pop(xml_out);  // Lsink_Lsrc
      pop(xml_out); //Pseudoscalars

    }

  // -------------------------



    QDPIO::cout << "Computed 8 basic pseudoscalars"  << endl;


  }



} // end of namespace
