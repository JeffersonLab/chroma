// $Id: quarkprop4_s.cc,v 3.7 2008-09-06 18:36:25 bjoo Exp $
/*! \file
 *  \brief Full quark propagator solver
 *
 *  Given a complete propagator as a source, this does all the inversions needed
 */

#include "chromabase.h"
#include "fermact.h"
#include "util/ferm/transf.h"
#include "actions/ferm/qprop/quarkprop4_s.h"
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
   * \param invParam inverter parameters ( Read )
   * \param RsdCG    CG (or MR) residual used here ( Read )
   * \param MaxCG    maximum number of CG iterations ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */

  template<typename T, typename P, typename Q>
  void quarkProp_a(LatticeStaggeredPropagator& q_sol, 
		   XMLWriter& xml_out,
		   const LatticeStaggeredPropagator& q_src,
		   const StaggeredTypeFermAct<T,P,Q>& S_f,
		   Handle< FermState<T,P,Q> > state,
		   const GroupXML_t& invParam,
		   QuarkSpinType quarkSpinType,
		   int& ncg_had)
  {
    START_CODE();

    push(xml_out, "QuarkProp4");

    ncg_had = 0;

    Handle<const SystemSolver<T> > qprop(S_f.qprop(state,invParam));

//  LatticeStaggeredFermion psi = zero;  // note this is ``zero'' and not 0

    // This version loops over all color and spin indices
    for(int color_source = 0; color_source < Nc; ++color_source)
    {
      QDPIO::cout<<"quarkprop_s:: doing color  : "<< color_source<<endl;

      LatticeStaggeredFermion psi = zero;  // note this is ``zero'' and not 0
      LatticeStaggeredFermion chi;

      // Extract a fermion source
      PropToFerm(q_src, chi, color_source);

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

      // Compute the propagator for given source color.
      {
	SystemSolverResults_t result = (*qprop)(psi,chi);
	ncg_had += result.n_count;

	push(xml_out,"Qprop");
	write(xml_out, "color_source", color_source);
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
      FermToProp(psi, q_sol, color_source);
    } /* end loop over color_source */

    pop(xml_out);

    END_CODE();
  }




  typedef LatticeStaggeredFermion LF;
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
  void quarkProp4(LatticeStaggeredPropagator& q_sol, 
		  XMLWriter& xml_out,
		  const LatticeStaggeredPropagator& q_src,	
		  const StaggeredTypeFermAct<LatticeStaggeredFermion, 
		  multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& S_f,
		  Handle< FermState<LatticeStaggeredFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> > > state,
		  const GroupXML_t& invParam,
		  QuarkSpinType quarkSpinType,
		  int& ncg_had)
  {
    quarkProp_a<LatticeStaggeredFermion,
      multi1d<LatticeColorMatrix>,
      multi1d<LatticeColorMatrix> >(q_sol, xml_out, q_src, S_f, state, invParam, 
				    quarkSpinType, ncg_had);
  }



  //! Given a complete propagator as a source, this does all the inversions needed
  /*! \ingroup qprop
   *
   * This routine is actually generic to all Staggered-like fermions
   *
   * \param q_sol    quark propagator ( Write )
   * \param q_src    source ( Read )
   * \param invParam inverter parameters ( Read )
   * \param ncg_had  number of CG iterations ( Write )
   */
  template<>
  void 
  StaggeredTypeFermAct<LatticeStaggeredFermion, 
		       multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >::quarkProp(
    LatticeStaggeredPropagator& q_sol, 
    XMLWriter& xml_out,
    const LatticeStaggeredPropagator& q_src,
    Handle< FermState<LatticeStaggeredFermion,
    multi1d<LatticeColorMatrix>,
    multi1d<LatticeColorMatrix> > > state,
    const GroupXML_t& invParam,
    QuarkSpinType quarkSpinType,
    int& ncg_had) const
  {
    quarkProp4(q_sol, xml_out, q_src, *this, state, invParam, quarkSpinType, ncg_had);
  }





  // Return a linear operator solver for this action to solve M*psi=chi 
  /*! \ingroup qprop */
  template<>
  LinOpSystemSolver<LF>*
  StaggeredTypeFermAct<LF,LCM,LCM>::invLinOp(Handle< FermState<LF,LCM,LCM> > state,
					     const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);
	
    // THIS NEEDS TO BE FIXED TO USE A PROPER MDAGM
    return TheLinOpStagFermSystemSolverFactory::Instance().createObject(invParam.id,
									paramtop,
									invParam.path,
									this->linOp(state));
  }


  //! Return a linear operator solver for this action to solve MdagM*psi=chi 
  /*! \ingroup qprop */
  template<>
  MdagMSystemSolver<LF>*
  StaggeredTypeFermAct<LF,LCM,LCM>::invMdagM(Handle< FermState<LF,LCM,LCM> > state,
					     const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    // THIS NEEDS TO BE FIXED TO USE A PROPER MDAGM
    return TheMdagMStagFermSystemSolverFactory::Instance().createObject(invParam.id,
									paramtop,
									invParam.path,
									this->linOp(state));
  }


  //! Return a linear operator solver for this action to solve (M+shift_i)*psi_i = chi 
  /*! \ingroup qprop */
  template<>
  LinOpMultiSystemSolver<LF>*
  StaggeredTypeFermAct<LF,LCM,LCM>::mInvLinOp(Handle< FermState<LF,LCM,LCM> > state,
					      const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    return TheLinOpStagFermMultiSystemSolverFactory::Instance().createObject(invParam.id,
									     paramtop,
									     invParam.path,
									     this->linOp(state));
  }


  //! Return a linear operator solver for this action to solve (MdagM+shift_i)*psi_i = chi 
  /*! \ingroup qprop */
  template<>
  MdagMMultiSystemSolver<LF>*
  StaggeredTypeFermAct<LF,LCM,LCM>::mInvMdagM(Handle< FermState<LF,LCM,LCM> > state,
					      const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    // THIS NEEDS TO BE FIXED TO USE A PROPER MDAGM
    return TheMdagMStagFermMultiSystemSolverFactory::Instance().createObject(invParam.id,
									 paramtop,
									 invParam.path,
									 this->linOp(state));
  }

  //! Return a linear operator solver for this action to solve (MdagM+shift_i)*psi_i = chi 
  /*! \ingroup qprop */
  template<>
  MdagMMultiSystemSolverAccumulate<LF>*
  StaggeredTypeFermAct<LF,LCM,LCM>::mInvMdagMAcc(Handle< FermState<LF,LCM,LCM> > state,
					      const GroupXML_t& invParam) const
  {
    std::istringstream  xml(invParam.xml);
    XMLReader  paramtop(xml);

    // THIS NEEDS TO BE FIXED TO USE A PROPER MDAGM
    return TheMdagMStagFermMultiSystemSolverAccumulateFactory::Instance().createObject(invParam.id,
									 paramtop,
									 invParam.path,
									 this->linOp(state));
  }


} // namespace Chroma
