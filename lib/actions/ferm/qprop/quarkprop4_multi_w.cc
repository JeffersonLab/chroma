// $Id: quarkprop4_multi_w.cc,v 1.4 2004-05-27 11:21:23 bjoo Exp $
/*! \file
 *  \brief Full quark propagator solver
 *
 *  Given a complete propagator as a source, this does all the inversions needed
 */

#include "chromabase.h"
#include "fermact.h"
#include "actions/ferm/qprop/quarkprop4_multi_w.h"
#include "util/ferm/transf.h"


template<typename C>
static 
void multiQuarkProp4_m(multi1d<LatticePropagator>& q_sol, 
		       XMLWriter& xml_out,
		       const LatticePropagator& q_src,
		       const C& S_f,
		       Handle<const ConnectState> state,
		       enum InvType invType,
		       const multi1d<Real>& masses,
		       const multi1d<Real>& RsdCG, 
		       int MaxCG, int& ncg_had)
{
  START_CODE("multiQuarkProp4");
  
  push(xml_out, "multiQuarkProp4");

  // Ensure that q_sol is adequate
  if ( q_sol.size() != masses.size() ) { 
    q_sol.resize(masses.size());
  }

  // Grab a RsdCG like array to pass down.
  // Set all elements to RsdCG
  ncg_had = 0;

  multi1d<LatticeFermion> psi(masses.size());

  // This version loops over all color and spin indices
  for(int color_source = 0; color_source < Nc; ++color_source)
  {
    for(int spin_source = 0; spin_source < Ns; ++spin_source)
    {
      LatticeFermion chi;

      // Extract a fermion source
      PropToFerm(q_src, chi, color_source, spin_source);

      // Use the last initial guess as the current initial guess

      /* 
       * Normalize the source in case it is really huge or small - 
       * a trick to avoid overflows or underflows
       */
      Real fact = Real(1) / sqrt(norm2(chi));
      chi *= fact;

      // Compute the propagator for given source color/spin.
      int n_count;

      // The psi-s are zeroed in multiQprop
      S_f.multiQprop(psi, masses, state, chi, invType, RsdCG, 
		     1, MaxCG, n_count);


      ncg_had += n_count;

      push(xml_out,"Qprop");
      write(xml_out, "RsdCG", RsdCG);
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

  END_CODE("multiQuarkProp4");
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

void multiQuarkProp4(multi1d<LatticePropagator>& q_sol, 
		XMLWriter& xml_out,
		const LatticePropagator& q_src,
		const OverlapFermActBase& S_f,
		Handle<const ConnectState> state,
		enum InvType invType,
		const multi1d<Real>& masses,
		const multi1d<Real>& RsdCG, 
		int MaxCG, int& ncg_had)
{
  multiQuarkProp4_m(q_sol, xml_out, q_src, S_f, state, invType, masses, 
	       RsdCG, MaxCG, ncg_had);
}
