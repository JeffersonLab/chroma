// $Id: quarkprop4_w.cc,v 1.20 2005-01-14 20:13:06 edwards Exp $
/*! \file
 *  \brief Full quark propagator solver
 *
 *  Given a complete propagator as a source, this does all the inversions needed
 */

#include "chromabase.h"
#include "fermact.h"
#include "util/ferm/transf.h"
#include "actions/ferm/qprop/quarkprop4_w.h"


namespace Chroma 
{
  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop
   *
   * This routine is actually generic to all Wilson-like fermions
   *
   * \param q_sol    quark propagator ( Write )
   * \param q_src    source ( Read )
   * \param invParam inverter parameters ( Read )
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */

  template<typename T>
  void quarkProp4_a(LatticePropagator& q_sol, 
		    XMLWriter& xml_out,
		    const LatticePropagator& q_src,
		    const FermionAction<T>& S_f,
		    Handle<const ConnectState> state,
		    const InvertParam_t& invParam,
		    bool nonRelProp,
		    int& ncg_had)
  {
    START_CODE();

    push(xml_out, "QuarkProp4");

    ncg_had = 0;

    Handle<const SystemSolver<T> > qprop(S_f.qprop(state,invParam));

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
//	S_f.qprop(psi, state, chi, invParam, n_count);
	int n_count = (*qprop)(psi,chi);
	ncg_had += n_count;

	push(xml_out,"Qprop");
	write(xml_out, "RsdCG", invParam.RsdCG);
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
   * \param invParam inverter parameters ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  void quarkProp4(LatticePropagator& q_sol, 
		  XMLWriter& xml_out,
		  const LatticePropagator& q_src,
		  const WilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >& S_f,
		  Handle<const ConnectState> state,
		  const InvertParam_t& invParam,
		  bool nonRelProp,
		  int& ncg_had)
  {
    quarkProp4_a<LatticeFermion>(q_sol, xml_out, q_src, S_f, state, invParam, nonRelProp, ncg_had);
  }


  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop
   *
   * This routine is actually generic to all Wilson-like fermions
   *
   * \param q_sol    quark propagator ( Write )
   * \param q_src    source ( Read )
   * \param invParam inverter parameters ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  void quarkProp4(LatticePropagator& q_sol, 
		  XMLWriter& xml_out,
		  const LatticePropagator& q_src,
		  const WilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >& S_f,
		  Handle<const ConnectState> state,
		  const InvertParam_t& invParam,
		  bool nonRelProp,
		  int& ncg_had)
  {
    quarkProp4_a<LatticeFermion>(q_sol, xml_out, q_src, S_f, state, invParam, nonRelProp, ncg_had);
  }


  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop
   *
   * This routine is actually generic to all Wilson-like fermions
   *
   * \param q_sol    quark propagator ( Write )
   * \param q_src    source ( Read )
   * \param invParam inverter parameters ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  template<>
  void 
  WilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >::quarkProp(LatticePropagator& q_sol, 
									      XMLWriter& xml_out,
									      const LatticePropagator& q_src,
									      Handle<const ConnectState> state,
									      const InvertParam_t& invParam,
									      bool nonRelProp,
									      int& ncg_had)
  {
    quarkProp4(q_sol, xml_out, q_src, *this, state, invParam, nonRelProp, ncg_had);
  }


  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop
   *
   * This routine is actually generic to all Wilson-like fermions
   *
   * \param q_sol    quark propagator ( Write )
   * \param q_src    source ( Read )
   * \param invParam inverter parameters ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  template<>
  void 
  WilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >::quarkProp(LatticePropagator& q_sol, 
										XMLWriter& xml_out,
										const LatticePropagator& q_src,
										Handle<const ConnectState> state,
										const InvertParam_t& invParam,
										bool nonRelProp,
										int& ncg_had)
  {
    quarkProp4(q_sol, xml_out, q_src, *this, state, invParam, nonRelProp, ncg_had);
  }

}; // namespace Chroma
