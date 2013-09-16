// $Id: quarkprop4_w.cc,v 3.8 2009-04-17 02:05:34 bjoo Exp $
/*! \file
 *  \brief Full quark propagator solver
 *
 *  Given a complete propagator as a source, this does all the inversions needed
 */

#include "wilstype_fermact_w.h"
#include "util/ferm/transf.h"
#include "actions/ferm/qprop/quarkprop4_w.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/multi_syssolver_linop_factory.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_factory.h"
#include "actions/ferm/invert/multi_syssolver_mdagm_accumulate_factory.h"


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
	{
	  SystemSolverResults_t result = (*qprop)(psi,chi);
	  ncg_had += result.n_count;

	  push(xml_out,"Qprop");
	  write(xml_out, "color_source", color_source);
	  write(xml_out, "spin_source", spin_source);
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
		  int& ncg_had)
  {
    quarkProp4_a<LF>(q_sol, xml_out, q_src, qprop, quarkSpinType, ncg_had);
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
    const GroupXML_t& invParam,
    QuarkSpinType quarkSpinType,
    int& ncg_had) const
  {
    Handle< SystemSolver<LF> > qprop(this->qprop(state,invParam));
    quarkProp4_a<LF>(q_sol, xml_out, q_src, qprop, quarkSpinType, ncg_had);
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
    const GroupXML_t& invParam,
    QuarkSpinType quarkSpinType,
    int& ncg_had) const
  {
    Handle< SystemSolver<LF> > qprop(this->qprop(state,invParam));
    quarkProp4_a<LF>(q_sol, xml_out, q_src, qprop, quarkSpinType, ncg_had);
  }



  //------------------------------------------------------------------------------------

  // Return a linear operator solver for this action to solve M*psi=chi 
  /*! \ingroup qprop */
  template<>
  LinOpSystemSolver<LF>*
  WilsonTypeFermAct<LF,LCM,LCM>::invLinOp(Handle< FermState<LF,LCM,LCM> > state,
					  const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);
	
    return TheLinOpFermSystemSolverFactory::Instance().createObject(invParam.id,
								    paramtop,
								    invParam.path,
								    state,
								    this->linOp(state));
  }


  //! Return a linear operator solver for this action to solve MdagM*psi=chi 
  /*! \ingroup qprop */
  template<>
  MdagMSystemSolver<LF>*
  WilsonTypeFermAct<LF,LCM,LCM>::invMdagM(Handle< FermState<LF,LCM,LCM> > state,
					  const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    return TheMdagMFermSystemSolverFactory::Instance().createObject(invParam.id,
								    paramtop,
								    invParam.path,
								    state,
								    this->linOp(state));
  }


  //! Return a linear operator solver for this action to solve (M+shift_i)*psi_i = chi 
  /*! \ingroup qprop */
  template<>
  LinOpMultiSystemSolver<LF>*
  WilsonTypeFermAct<LF,LCM,LCM>::mInvLinOp(Handle< FermState<LF,LCM,LCM> > state,
					   const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    return TheLinOpFermMultiSystemSolverFactory::Instance().createObject(invParam.id,
									 paramtop,
									 invParam.path,
									 this->linOp(state));
  }


  //! Return a linear operator solver for this action to solve (MdagM+shift_i)*psi_i = chi 
  /*! \ingroup qprop */
  template<>
  MdagMMultiSystemSolver<LF>*
  WilsonTypeFermAct<LF,LCM,LCM>::mInvMdagM(Handle< FermState<LF,LCM,LCM> > state,
					   const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    return TheMdagMFermMultiSystemSolverFactory::Instance().createObject(invParam.id,
									 paramtop,
									 invParam.path,
									 state,
									 this->linOp(state));
  }

  //! Return a linear operator solver for this action to solve (MdagM+shift_i)*psi_i = chi 
  /*! \ingroup qprop */
  template<>
  MdagMMultiSystemSolverAccumulate<LF>*
  WilsonTypeFermAct<LF,LCM,LCM>::mInvMdagMAcc(Handle< FermState<LF,LCM,LCM> > state,
					   const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    return TheMdagMFermMultiSystemSolverAccumulateFactory::Instance().createObject(invParam.id,
									 paramtop,
									 invParam.path,
									 this->linOp(state));
  }



  //------------------------------------------------------------------------------------

  // Return a linear operator solver for this action to solve M*psi=chi 
  /*! \ingroup qprop */
  template<>
  LinOpSystemSolverArray<LF>*
  WilsonTypeFermAct5D<LF,LCM,LCM>::invLinOp(Handle< FermState<LF,LCM,LCM> > state,
					    const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);
	
    return TheLinOpFermSystemSolverArrayFactory::Instance().createObject(invParam.id,
									 paramtop,
									 invParam.path,
									 this->linOp(state));
  }


  //! Return a linear operator solver for this action to solve MdagM*psi=chi 
  /*! \ingroup qprop */
  template<>
  MdagMSystemSolverArray<LF>*
  WilsonTypeFermAct5D<LF,LCM,LCM>::invMdagM(Handle< FermState<LF,LCM,LCM> > state,
					    const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    return TheMdagMFermSystemSolverArrayFactory::Instance().createObject(invParam.id,
									 paramtop,
									 invParam.path,
									 this->linOp(state));
  }



  // Return a linear operator solver for this action to solve M*psi=chi 
  /*! \ingroup qprop */
  template<>
  LinOpSystemSolverArray<LF>*
  WilsonTypeFermAct5D<LF,LCM,LCM>::invLinOpPV(Handle< FermState<LF,LCM,LCM> > state,
					      const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);
	
    return TheLinOpFermSystemSolverArrayFactory::Instance().createObject(invParam.id,
									 paramtop,
									 invParam.path,
									 this->linOpPV(state));
  }


  //! Return a linear operator solver for this action to solve PV^dag*PV*psi=chi 
  /*! \ingroup qprop */
  template<>
  MdagMSystemSolverArray<LF>*
  WilsonTypeFermAct5D<LF,LCM,LCM>::invMdagMPV(Handle< FermState<LF,LCM,LCM> > state,
					      const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    return TheMdagMFermSystemSolverArrayFactory::Instance().createObject(invParam.id,
									 paramtop,
									 invParam.path,
									 this->linOpPV(state));
  }


  //! Return a linear operator solver for this action to solve (MdagM+shift_i)*psi_i = chi 
  /*! \ingroup qprop */
  template<>
  MdagMMultiSystemSolverArray<LF>*
  WilsonTypeFermAct5D<LF,LCM,LCM>::mInvMdagM(Handle< FermState<LF,LCM,LCM> > state,
					     const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    return TheMdagMFermMultiSystemSolverArrayFactory::Instance().createObject(invParam.id,
									      paramtop,
									      invParam.path,
									      state,
									      lMdagM(state));
  }

  //! Return a linear operator solver for this action to solve (MdagM+shift_i)*psi_i = chi 
  /*! \ingroup qprop */
  template<>
  MdagMMultiSystemSolverAccumulateArray<LF>*
  WilsonTypeFermAct5D<LF,LCM,LCM>::mInvMdagMAcc(Handle< FermState<LF,LCM,LCM> > state,
					     const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    return TheMdagMFermMultiSystemSolverAccumulateArrayFactory::Instance().createObject(invParam.id,
									      paramtop,
									      invParam.path,
									      lMdagM(state));
  }


  //! Return a linear operator solver for this action to solve (MdagM+shift_i)*psi_i = chi 
  /*! \ingroup qprop */
  template<>
  MdagMMultiSystemSolverArray<LF>*
  WilsonTypeFermAct5D<LF,LCM,LCM>::mInvMdagMPV(Handle< FermState<LF,LCM,LCM> > state,
					       const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    Handle< LinearOperatorArray<LF> > PV(this->linOpPV(state));
    Handle< LinearOperatorArray<LF> > MdagM(new MdagMLinOpArray<LF>(PV));

    return TheMdagMFermMultiSystemSolverArrayFactory::Instance().createObject(
      invParam.id,
      paramtop,
      invParam.path,
      state,
      MdagM);
  }

  //! Return a linear operator solver for this action to solve (MdagM+shift_i)*psi_i = chi 
  /*! \ingroup qprop */
  template<>
  MdagMMultiSystemSolverAccumulateArray<LF>*
  WilsonTypeFermAct5D<LF,LCM,LCM>::mInvMdagMPVAcc(Handle< FermState<LF,LCM,LCM> > state,
					       const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    Handle< LinearOperatorArray<LF> > PV(this->linOpPV(state));
    Handle< LinearOperatorArray<LF> > MdagM(new MdagMLinOpArray<LF>(PV));

    return TheMdagMFermMultiSystemSolverAccumulateArrayFactory::Instance().createObject(
      invParam.id,
      paramtop,
      invParam.path,
      MdagM);
  }


} // namespace Chroma
