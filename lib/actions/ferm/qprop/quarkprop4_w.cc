// $Id: quarkprop4_w.cc,v 3.1 2006-06-11 06:30:33 edwards Exp $
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
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */

  template<typename T>
  void quarkProp4_a(LatticePropagator& q_sol, 
		    XMLWriter& xml_out,
		    const LatticePropagator& q_src,
		    Handle< SystemSolver<T> > qprop,
		    QuarkSpinType quarkSpinType,
		    int numRetries,
		    int& ncg_had)
  {
    START_CODE();

    QDPIO::cout << "Entering quarkProp4" << endl;
    push(xml_out, "QuarkProp4");

    ncg_had = 0;

    int start_spin;
    int end_spin;

    switch (quarkSpinType)
    {
    case QUARK_SPIN_TYPE_FULL:
      start_spin = 0;
      end_spin = Ns;
      break;

    case QUARK_SPIN_TYPE_UPPER:
      start_spin = 0;
      end_spin = Ns/2;
      break;

    case QUARK_SPIN_TYPE_LOWER:
      start_spin = Ns/2;
      end_spin = Ns;
      break;
    }

//  LatticeFermion psi = zero;  // note this is ``zero'' and not 0

    // This version loops over all color and spin indices
    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      for(int spin_source = start_spin; spin_source < end_spin; ++spin_source)
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
	Real fact = 1.0;
	Real nrm = sqrt(norm2(chi));
	if (toFloat(nrm) != 0.0)
	  fact /= nrm;

	// Rescale
	chi *= fact;

	// Compute the propagator for given source color/spin.
	for(int ntry=0; ntry < numRetries; ++ntry)
	{
	  SystemSolverResults_t result = (*qprop)(psi,chi);
	  ncg_had += result.n_count;

	  push(xml_out,"Qprop");
	  write(xml_out, "color_source", color_source);
	  write(xml_out, "spin_source", spin_source);
	  write(xml_out, "ntry", ntry);
	  write(xml_out, "n_count", result.n_count);
	  write(xml_out, "resid", result.resid);
	  pop(xml_out);
	}

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


    switch (quarkSpinType)
    {
    case QUARK_SPIN_TYPE_FULL:
      // Do nothing here
      break;

    case QUARK_SPIN_TYPE_UPPER:
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
	for(int spin_source = Ns/2; spin_source < Ns; ++spin_source) 
	{ 
	  int copyfrom = spin_source - Ns/2;
	  LatticeFermion psi;

	  PropToFerm(q_sol, psi, color_source, copyfrom);
	  FermToProp(LatticeFermion(-psi), q_sol, color_source, spin_source);
	}
      }
    }
    break;

    case QUARK_SPIN_TYPE_LOWER:
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
	for(int spin_source = 0; spin_source < Ns/2; ++spin_source) 
	{ 
	  int copyfrom = spin_source + Ns/2;
	  LatticeFermion psi;

	  PropToFerm(q_sol, psi, color_source, copyfrom);
	  // There is no need for (-) in the lower component case (KNO)
	  FermToProp(LatticeFermion(psi), q_sol, color_source, spin_source);
	}
      }
    }
    break;
    }  // end switch(quarkSpinType)

    pop(xml_out);
    QDPIO::cout << "Exiting quarkProp4" << endl;

    END_CODE();
  }


  typedef LatticeFermion LF;
  typedef multi1d<LatticeColorMatrix> LCM;


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
		  Handle< SystemSolver<LF> > qprop,
		  QuarkSpinType quarkSpinType,
		  int numRetries,
		  int& ncg_had)
  {
    quarkProp4_a<LF>(q_sol, xml_out, q_src, qprop, quarkSpinType, numRetries, ncg_had);
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
  WilsonTypeFermAct<LF,LCM,LCM>::quarkProp(
    LatticePropagator& q_sol, 
    XMLWriter& xml_out,
    const LatticePropagator& q_src,
    Handle< FermState<LF,LCM,LCM> > state,
    const InvertParam_t& invParam,
    QuarkSpinType quarkSpinType,
    int numRetries,
    int& ncg_had) const
  {
    Handle< SystemSolver<LF> > qprop(this->qprop(state,invParam));
    quarkProp4_a<LF>(q_sol, xml_out, q_src, qprop, quarkSpinType, numRetries, ncg_had);
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
  WilsonTypeFermAct5D<LF,LCM,LCM>::quarkProp(
    LatticePropagator& q_sol, 
    XMLWriter& xml_out,
    const LatticePropagator& q_src,
    Handle< FermState<LF,LCM,LCM> > state,
    const InvertParam_t& invParam,
    QuarkSpinType quarkSpinType,
    int numRetries,
    int& ncg_had) const
  {
    Handle< SystemSolver<LF> > qprop(this->qprop(state,invParam));
    quarkProp4_a<LF>(q_sol, xml_out, q_src, qprop, quarkSpinType, numRetries, ncg_had);
  }

} // namespace Chroma
