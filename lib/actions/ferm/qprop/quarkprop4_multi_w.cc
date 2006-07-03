// $Id: quarkprop4_multi_w.cc,v 3.2 2006-07-03 15:26:09 edwards Exp $
/*! \file
 *  \brief Full quark propagator solver
 *
 *  Given a complete propagator as a source, this does all the inversions needed
 */

#include "chromabase.h"
#include "fermact.h"
#include "util/ferm/transf.h"
#include "actions/ferm/fermacts/overlap_fermact_base_w.h"
#include "io/io.h"

namespace Chroma { 
  template<class T, class C>
  static 
  void multiQuarkProp4_m(multi1d<LatticePropagator>& q_sol, 
			 XMLWriter& xml_out,
			 const LatticePropagator& q_src,
			 const C& S_f,
			 Handle< FermState<LatticeFermion,
			 multi1d<LatticeColorMatrix>,
			 multi1d<LatticeColorMatrix> > > state,
			 const multi1d<Real>& masses,
			 const GroupXML_t& invParam,
			 int n_soln,
			 int& ncg_had)
  {
    START_CODE();
    
    push(xml_out, "multiQuarkProp4");
    
    // Ensure that q_sol is adequate
    if ( q_sol.size() != masses.size() ) { 
      q_sol.resize(masses.size());
    }
    
    // Grab a RsdCG like array to pass down.
    // Set all elements to RsdCG
    ncg_had = 0;
    
    multi1d<T> psi(masses.size());
    
    // This version loops over all color and spin indices
    for(int color_source = 0; color_source < Nc; ++color_source)
      {
	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	  {
	    T chi;
	    
	    // Extract a fermion source
	    PropToFerm(q_src, chi, color_source, spin_source);
	    
	    // Use the last initial guess as the current initial guess
	    
	    /* 
	     * Normalize the source in case it is really huge or small - 
	     * a trick to avoid overflows or underflows
	     */
	    Real fact = 1.0;
	    Real nrm = sqrt(norm2(chi));
	    if (toFloat(nrm) != 0.0)
	      fact /= nrm;

	    // Rescale
	    chi *= fact;
	    
	    // Compute the propagator for given source color/spin.
	    int n_count;
	    
	    // The psi-s are zeroed in multiQprop
	    S_f.multiQprop(psi, masses, state, chi, invParam, n_soln, n_count);
	    
	    ncg_had += n_count;
	    
	    push(xml_out,"Qprop");
//	    write(xml_out, "RsdCG", invParam.RsdCG);
	    write(xml_out, "n_count", n_count);
	    pop(xml_out);
	    
	    // Unnormalize the source following the inverse of the normalization above
	    fact = Real(1) / fact;
	    
	    for(int i=0; i < masses.size(); i++) { 
	      psi[i] *= fact;
	      /*
	       * Move the solution to the appropriate components
	       * of quark propagator.
	       */
	      FermToProp(psi[i], q_sol[i], color_source, spin_source);
	    }
	    
	  }	/* end loop over spin_source */
      } /* end loop over color_source */
    
    pop(xml_out);
    
    END_CODE();
  }
  
//! Given a complete propagator as a source, this does all the inversions needed
/*! \ingroup qprop
 *
 * This routine is actually generic to all Wilson-like fermions
 *
 * \param q_sol    quark propagator ( Write )
 * \param q_src    source ( Read )
 * \param invType  inverter type ( Read (
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

  void OverlapFermActBase::multiQuarkProp(multi1d<LatticePropagator>& q_sol, 
					  XMLWriter& xml_out,
					  const LatticePropagator& q_src,
					  Handle< FermState<LatticeFermion,
					  multi1d<LatticeColorMatrix>,
					  multi1d<LatticeColorMatrix> > > state,
					  const multi1d<Real>& masses,
					  const GroupXML_t& invParam,
					  const int n_soln,
					  int& ncg_had)
  {
    multiQuarkProp4_m<LatticeFermion, OverlapFermActBase>(q_sol, 
							  xml_out, 
							  q_src, 
							  *this, 
							  state, 
							  masses, 
							  invParam, 
							  n_soln,
							  ncg_had);
  }

};
