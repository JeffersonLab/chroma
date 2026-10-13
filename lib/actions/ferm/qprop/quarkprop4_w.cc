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

    QDPIO::cout << "Entering quarkProp4 - MRHS interface" << std::endl;
    push(xml_out, "QuarkProp4");

    ncg_had = 0;

    int start_spin;
    int end_spin;
		int num_spin;

    switch (quarkSpinType)
    {
    case QUARK_SPIN_TYPE_FULL:
      start_spin = 0;
      end_spin = Ns;
			num_spin = Ns;
      break;

    case QUARK_SPIN_TYPE_UPPER:
      start_spin = 0;
      end_spin = Ns/2;
			num_spin = Ns/2;
      break;

    case QUARK_SPIN_TYPE_LOWER:
      start_spin = Ns/2;
      end_spin = Ns;
			num_spin = Ns/2;
      break;
    }


		{ 
			multi1d<Double> norm_chi(Nc*num_spin);
			multi1d<Double> fact(Nc*num_spin);
			std::vector< std::shared_ptr<const LatticeFermion> > chi_ptrs(Nc*num_spin);
			std::vector< std::shared_ptr<LatticeFermion> > psi_ptrs(Nc*num_spin);
			// This version loops over all color and spin indices
			int idx=0;
			for(int color_source = 0; color_source < Nc; ++color_source)
			{
				for(int spin_source = start_spin; spin_source < end_spin; ++spin_source)
				{
					psi_ptrs[idx] = std::make_shared<LatticeFermion>(zero);
					

					// Extract a fermion source
					// Due to the vaguaries of initializing a std::shared<const T>
					// We go via a temporary.
					LatticeFermion tmp;
					PropToFerm(q_src, tmp, color_source, spin_source);

					// Normalize temporary 
					norm_chi[idx] = sqrt(norm2(tmp));
					fact[idx] = toDouble(1)/norm_chi[idx];
					tmp *= fact[idx];
				
					// Create the RHS 	
					chi_ptrs[idx] = std::make_shared<const LatticeFermion>(tmp);

					// Update Index
					idx++;
				}
			}

			// Do the MultiRHS solve
			//
			// Convention: In true multiRHS solve only solution 0 will have non-zero
			// n-count for now. That way accumulating ncg_had by adding 0s potentially
			// will work.
			std::vector<SystemSolverResults_t> results = (*qprop)(psi_ptrs, chi_ptrs);

			// Accumulate ncg_had and restore solution into solution prop	
			ncg_had = 0;
			for(int idx=0; idx < Nc*num_spin; idx++) {

				// Undo rescale by multiplying by 1/fact = norm_chi[idx]
				*(psi_ptrs[idx]) *= norm_chi[idx]; 

				// break colorspin index into color and spin indices. 
				int spin_idx = idx%num_spin + start_spin;
				int col_idx =idx/num_spin; 

				// Insert  solution into propagator
				FermToProp(*(psi_ptrs[idx]), q_sol, col_idx, spin_idx);

				// Accumulate ncg_had. This will be correct if we follow
				// the convention that true mrhs solvers return only a count
				// in results[0].n_count and keep all others as zero
				// Fake MRHS solvers (which loop over sources) can fill out 
				// an accurate iteration count for each solve. 
				ncg_had += results[idx].n_count;
				push(xml_out,"Qprop");
				write(xml_out, "color_source", col_idx);
				write(xml_out, "spin_source", spin_idx);
				write(xml_out, "n_count", results[idx].n_count);
				write(xml_out, "resid", results[idx].resid);
				pop(xml_out);

			} /* end loop over solutions */
		} // psis, chis etc go away here. 

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
				for(int color_source = 0; color_source < Nc ; ++color_source) {
					for(int spin_source = Ns/2; spin_source < Ns; ++spin_source) { int copyfrom = spin_source - Ns/2;
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
				for(int color_source = 0; color_source < Nc ; ++color_source) {
					for(int spin_source = 0; spin_source < Ns/2; ++spin_source) { 
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
    QDPIO::cout << "Exiting quarkProp4" << std::endl;

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
    QDPIO::cout << "In quarkProp()" << std::endl;
    StopWatch swatch;
    swatch.start();
    Handle< SystemSolver<LF> > qprop(this->qprop(state,invParam));
    swatch.stop();
    QDPIO::cout << "Creating qprop took " << swatch.getTimeInSeconds() 
		<< "sec " << std::endl;
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
									 state,
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
									 state,
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
									 state,
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
									 state,
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
