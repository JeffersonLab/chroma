// -*- C++ -*-
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2 using CG
 */

#ifndef __multi_syssolver_mdagm_cg_chrono_h__
#define __multi_syssolver_mdagm_cg_chrono_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/multi_syssolver_mdagm.h"
#include "actions/ferm/invert/multi_syssolver_cg_params.h"
#include "actions/ferm/invert/invcg2.h"
#include "actions/ferm/linop/lopishift.h"
#include "update/molecdyn/predictor/mre_extrap_predictor.h"
#include "init/chroma_init.h"

namespace Chroma
{

  //! CG2 system solver namespace
  namespace MdagMMultiSysSolverCGChronoEnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a CG2 system. Here, the operator is NOT assumed to be hermitian
  /*! \ingroup invert
   */
  template<typename T>
  class MdagMMultiSysSolverCGChrono : public MdagMMultiSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    MdagMMultiSysSolverCGChrono(Handle< LinearOperator<T> > A_,
			  const MultiSysSolverCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~MdagMMultiSysSolverCGChrono() {}

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (multi1d<T>& psi, const multi1d<Real>& shifts, const T& chi) const
      {
	START_CODE();

	multi1d<Real> RsdCG(shifts.size());
	if (invParam.RsdCG.size() == 1)
	{
	  RsdCG = invParam.RsdCG[0];
	}
	else if (invParam.RsdCG.size() == RsdCG.size())
	{
	  RsdCG = invParam.RsdCG;
	}
	else
	{
	  QDPIO::cerr << "MdagMMultiSysSolverCGChrono: shifts incompatible" << endl;
	  QDP_abort(1);
	}
	MinimalResidualExtrapolation4DChronoPredictor<T> chrono(shifts.size());
	SystemSolverResults_t res;
	for(int i=shifts.size(); i >=0; i--) { 
	  Real rshift = sqrt(shifts[i]);
	  Handle<LinearOperator<T> > Ms(new lopishift<T,Real>(A, rshift));
	  chrono.predictX(psi[i], *Ms, chi);
	  res = InvCG2(*Ms, chi, psi[i], RsdCG[i], invParam.MaxCG);
	  chrono.newXVector(psi[i]);
	}
#if 1
	XMLFileWriter& log = Chroma::getXMLLogInstance();
	push(log, "MultiCG");
	write(log, "shifts", shifts);
	write(log, "RsdCG", RsdCG);
	write(log, "n_count", res.n_count);
	Double chinorm=norm2(chi, A->subset());
	multi1d<Double> r_rel(shifts.size());

	for(int i=0; i < shifts.size(); i++) { 
	   T tmp1,tmp2;
	   (*A)(tmp1, psi[i], PLUS);
	   (*A)(tmp2, tmp1, MINUS);  // tmp2 = A^\dagger A psi
	   tmp2[ A->subset() ] +=  shifts[i]* psi[i]; // tmp2 = ( A^\dagger A + shift_i ) psi
	   T r;
	   r[ A->subset() ] = chi - tmp2;
	   r_rel[i] = sqrt(norm2(r, A->subset())/chinorm );
	}
	write(log, "ResidRel", r_rel);
	pop(log);
#endif
	END_CODE();

	return res;
      }


  private:
    // Hide default constructor
    MdagMMultiSysSolverCGChrono() {}

    Handle< LinearOperator<T> > A;
    MultiSysSolverCGParams invParam;
  };


} // End namespace

#endif 

