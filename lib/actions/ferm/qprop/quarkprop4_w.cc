// $Id: quarkprop4_w.cc,v 1.15 2004-07-28 02:38:02 edwards Exp $
/*! \file
 *  \brief Full quark propagator solver
 *
 *  Given a complete propagator as a source, this does all the inversions needed
 */

#include "chromabase.h"
#include "fermact.h"
#include "actions/ferm/qprop/quarkprop4_w.h"
#include "util/ferm/transf.h"

using namespace QDP;

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

template<typename T, template<class> class C>
static 
void quarkProp4_a(LatticePropagator& q_sol, 
		  XMLWriter& xml_out,
		  const LatticePropagator& q_src,
		  const C<T>& S_f,
		  Handle<const ConnectState> state,
		  enum InvType invType,
		  const Real& RsdCG, 
		  int MaxCG, 
		  bool nonRelProp,
		  int& ncg_had)
{
  START_CODE();

  push(xml_out, "QuarkProp4");

  ncg_had = 0;

  int max_spin = (nonRelProp) ? (Ns/2) : Ns;

//  LatticeFermion psi = zero;  // note this is ``zero'' and not 0

  // This version loops over all color and spin indices
  for(int color_source = 0; color_source < Nc; ++color_source)
  {
    for(int spin_source = 0; spin_source < max_spin; ++spin_source)
    {
      LatticeFermion psi = zero;  // note this is ``zero'' and not 0
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

      S_f.qprop(psi, state, chi, invType, RsdCG, MaxCG, n_count);
      ncg_had += n_count;

      push(xml_out,"Qprop");
      write(xml_out, "RsdCG", RsdCG);
      write(xml_out, "n_count", n_count);
      pop(xml_out);

      // Unnormalize the source following the inverse of the normalization above
      fact = Real(1) / fact;
      psi *= fact;

      /*
       * Move the solution to the appropriate components
       * of quark propagator.
       */
      FermToProp(psi, q_sol, color_source, spin_source);
    }	/* end loop over spin_source */
  } /* end loop over color_source */


  if ( nonRelProp )
  {
    /* Since this is a non-relativistic prop 
     * negate the quark props 'lower' components
     * This is because I should have only done a half inversion 
     * on non relativistic channels, where the last two columns of the 
     * source MUST be the negation of the first two columns. 
     * Hence the last two columns of the solution must also be 
     * negations of the first two columns. The half inversion itself
     * has not put in the minus sign, it just copied the columns.
     * The post multiply by Gamma_5 adds in the required - sign 
     * in the last two columns 
     */ 
    /* Apply Gamma_5 = Gamma(15) by negating the fermion extracted */
    for(int color_source = 0; color_source < Nc ; ++color_source) 
    {
      for(int spin_source = max_spin; spin_source < Ns; ++spin_source) 
      { 
	int copyfrom = spin_source - max_spin;
	LatticeFermion psi;

	PropToFerm(q_sol, psi, color_source, copyfrom);
	FermToProp(LatticeFermion(-psi), q_sol, color_source, spin_source);
      }
    }
  }

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

void quarkProp4(LatticePropagator& q_sol, 
		XMLWriter& xml_out,
		const LatticePropagator& q_src,
		const WilsonTypeFermAct<LatticeFermion>& S_f,
		Handle<const ConnectState> state,
		enum InvType invType,
		const Real& RsdCG, 
		int MaxCG, 
		bool nonRelProp,
		int& ncg_had)
{
  quarkProp4_a(q_sol, xml_out, q_src, S_f, state, invType, RsdCG, MaxCG, nonRelProp, ncg_had);
}


//! Preconditioned version for overlap
void quarkProp4(LatticePropagator& q_sol, 
		XMLWriter& xml_out,
		const LatticePropagator& q_src,
		const OverlapFermActBase& S_f,
		Handle<const ConnectState> state,
		enum InvType invType,
		const Real& RsdCG, 
		const Real& RsdCGPrec,
		int MaxCG, 
		int MaxCGPrec,
		bool nonRelProp,
		int& ncg_had)
{
  START_CODE();

  push(xml_out, "QuarkProp4");

  ncg_had = 0;

  int max_spin = (nonRelProp) ? (Ns/2) : Ns;

//  LatticeFermion psi = zero;  // note this is ``zero'' and not 0

  // This version loops over all color and spin indices
  for(int color_source = 0; color_source < Nc; ++color_source)
  {
    for(int spin_source = 0; spin_source < max_spin; ++spin_source)
    {
      LatticeFermion psi = zero;  // note this is ``zero'' and not 0
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

      S_f.qprop(psi, state, chi, invType, RsdCG, RsdCGPrec, MaxCG, MaxCGPrec, n_count);
      ncg_had += n_count;

      push(xml_out,"Qprop");
      write(xml_out, "RsdCG", RsdCG);
      write(xml_out, "n_count", n_count);
      pop(xml_out);

      // Unnormalize the source following the inverse of the normalization above
      fact = Real(1) / fact;
      psi *= fact;

      /*
       * Move the solution to the appropriate components
       * of quark propagator.
       */
      FermToProp(psi, q_sol, color_source, spin_source);
    }	/* end loop over spin_source */
  } /* end loop over color_source */


  if ( nonRelProp )
  {
    /* Since this is a non-relativistic prop 
     * negate the quark props 'lower' components
     * This is because I should have only done a half inversion 
     * on non relativistic channels, where the last two columns of the 
     * source MUST be the negation of the first two columns. 
     * Hence the last two columns of the solution must also be 
     * negations of the first two columns. The half inversion itself
     * has not put in the minus sign, it just copied the columns.
     * The post multiply by Gamma_5 adds in the required - sign 
     * in the last two columns 
     */ 
    /* Apply Gamma_5 = Gamma(15) by negating the fermion extracted */
    for(int color_source = 0; color_source < Nc ; ++color_source) 
    {
      for(int spin_source = max_spin; spin_source < Ns; ++spin_source) 
      { 
	int copyfrom = spin_source - max_spin;
	LatticeFermion psi;

	PropToFerm(q_sol, psi, color_source, copyfrom);
	FermToProp(LatticeFermion(-psi), q_sol, color_source, spin_source);
      }
    }
  }

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

void quarkProp4(LatticePropagator& q_sol, 
		XMLWriter& xml_out,
		const LatticePropagator& q_src,
		const WilsonTypeFermAct< multi1d<LatticeFermion> >& S_f,
		Handle<const ConnectState> state,
		enum InvType invType,
		const Real& RsdCG, 
		int MaxCG, 
		bool nonRelProp,
		int& ncg_had)
{
  quarkProp4_a(q_sol, xml_out, q_src, S_f, state, invType, RsdCG, MaxCG, nonRelProp, ncg_had);
}


#include "actions/ferm/fermacts/prec_dwf_fermact_base_array_w.h"

//! Given a complete propagator as a source, this does all the inversions needed
/*! \ingroup qprop
 *
 * This routine is actually generic to Domain Wall fermions (Array) fermions
 *
 * \param q_sol    quark propagator ( Write )
 * \param q_src    source ( Read )
 * \param invType  inverter type ( Read (
 * \param RsdCG    CG (or MR) residual used here ( Read )
 * \param MaxCG    maximum number of CG iterations ( Read )
 * \param ncg_had  number of CG iterations ( Write )
 */

void quarkProp4(LatticePropagator& q_sol, 
		XMLWriter& xml_out,
		const LatticePropagator& q_src,
		const EvenOddPrecDWFermActBaseArray<LatticeFermion>& S_f,
		Handle<const ConnectState> state,
		enum InvType invType,
		const Real& RsdCG, 
		int MaxCG, 
		bool nonRelProp,
		int& ncg_had)
{
  quarkProp4_a(q_sol, xml_out, q_src, S_f, state, invType, RsdCG, MaxCG, nonRelProp, ncg_had);
}

