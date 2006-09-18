// -*- C++ -*-
// $Id: prec_dwf_qprop_array_cg_dwf_w.h,v 3.7 2006-09-18 18:23:35 bjoo Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#ifndef __prec_dwf_qprop_array_cg_dwf_w_h__
#define __prec_dwf_qprop_array_cg_dwf_w_h__

#include "fermact.h"
#include "io/aniso_io.h"
#include "actions/ferm/invert/syssolver_cg_params.h"

namespace Chroma
{

  //! SSE Propagator DWF qpropT
  /*! \ingroup qprop
   *
   * Propagator solver for DWF fermions
   */
  template< typename SinglePrecSolver, typename DoublePrecSolver >
  class CGDWFQpropT : public SystemSolverArray<LatticeFermion>
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;


    //! Alternative constructor for compatibility
    /*!
     * \param m_q_       quark mass ( Read )
     */
    CGDWFQpropT(Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A_,
		 Handle< LinOpSystemSolverArray<T> > invA_,   // throw away
		 Handle< FermState<T,P,Q> > state_, 
		 const Real& OverMass_,
		 const Real& Mass_,
		 const AnisoParam_t& anisoParam_,
		 const GroupXML_t& invParam_) : 
      A(A_), OverMass(OverMass_), Mass(Mass_), 
      N5(A->size()), anisoParam(anisoParam_)
      {init(state_, invParam_);}

    //! Need a real destructor
    ~CGDWFQpropT() {fini();}

    //! Expected length of array index
    int size() const {return N5;}

    //! Return the subset on which the operator acts
    const OrderedSubset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    SystemSolverResults_t operator() (multi1d<LatticeFermion>& psi, const multi1d<LatticeFermion>& chi) const
    {
      QDPIO::cout << "entering CGDWFQpropT::operator()" << endl;
      
      START_CODE();



      SystemSolverResults_t res;
      
      //    init();   // only needed because 2 qpropT might be active - SSE CG does not allow this
      
      Real ff = where(anisoParam.anisoP, anisoParam.nu / anisoParam.xi_0, Real(1));
      
      // Apply SSE inverter
      Real   a5  = 1;
      double M5  = toDouble(1 + a5*(1 + (Nd-1)*ff - OverMass));
      double m_f = toDouble(Mass);
      double rsd = toDouble(invParam.RsdCG);
      double rsd_sq = rsd * rsd;
      int    max_iter = invParam.MaxCG;
      double out_eps;

     
      // Do single prec solve first
      {
	QDPIO::cout << "CGDWFQpropT: Doing single prec solution" << endl;
	SinglePrecSolver single_prec_solver(lattice_size.slice(),
					    &u,
					    &v);
	
	single_prec_solver.cgSolver(psi, M5, m_f, 
				    chi, psi, rsd_sq, max_iter, out_eps, res.n_count);
      }


      // Do double prec solve second
      {
	QDPIO::cout << "CGDWFQpropT: Doing double prec restart" << endl;
	DoublePrecSolver double_prec_solver(lattice_size.slice(), 
					    &u,
					    &v);

	int n_count2;	
	double_prec_solver.cgSolver(psi, M5, m_f, 
				    chi, psi, rsd_sq, max_iter, out_eps, n_count2);

	res.n_count += n_count2;
      }

      // Compute residual
      {
	multi1d<LatticeFermion>  r(N5);
	A->unprecLinOp(r, psi, PLUS);
	r -= chi;
	res.resid = sqrt(norm2(r));
      }
      
      QDPIO::cout << "exiting CGDWFQpropT::operator()" << endl;
      
      END_CODE();
      
      return res;
    }

  protected:
    //! Private internal initializer: 
    //  Parse Params, 
    //  Prepare U field and V field
    //  Prepare lattice size
    void init(Handle< FermState<T,P,Q> > state, const GroupXML_t& inv)
    {

      // Parse params an
      QDPIO::cout << "entering CGDWFQpropT::init" << endl;
      
      if (Nd != 4 || Nc != 3) {
	
	QDPIO::cerr << "CGDWFQpropT: only supports Nd=4 and Nc=3" << endl;
	QDP_abort(1);
	
      }

      // Read the XML for the CG params
      try {
	
	std::istringstream  is(inv.xml);
	XMLReader  paramtop(is);
	
	read(paramtop, inv.path, invParam);
      }
      catch (const std::string& e) {
	
	QDPIO::cerr << "CGDWFQpropT: only support a CG inverter" << endl;
	QDP_abort(1);
      }


      lattice_size.resize(Nd+1);
      lattice_size[Nd] = N5;
      for(int i=0; i < Nd; ++i)
	lattice_size[i] = Layout::lattSize()[i];
      
      if (N5 % 4 != 0) {
	
	QDPIO::cerr << "SSE qpropT only N5 that is multiple of 2" << endl;
	QDP_abort(1);
      }

      // Transform the U fields
      u.resize(Nd);
      v.resize(Nd);

      u = state->getLinks();
      
      if (anisoParam.anisoP) {
	
	Real ff = where(anisoParam.anisoP, anisoParam.nu / anisoParam.xi_0, Real(1));
	
	// Rescale the u fields by the anisotropy
	for(int mu=0; mu < u.size(); ++mu) {
	  
	  if (mu != anisoParam.t_dir)
	    u[mu] *= ff;
	}
      }
      
      // Construct shifted gauge field
      for (int i = 0; i < Nd; i++)
	v[i] = shift(u[i], -1, i); // as viewed from the destination
      QDPIO::cout << "exiting CGDWFQpropT::init" << endl;
    }

    //! Private internal destructor
    void fini() {
      QDPIO::cout << "CGDWFQpropT: calling destructor" << endl;
    }
      
  private:
    Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A;
    Real OverMass;
    Real Mass;
    int  N5;
    AnisoParam_t anisoParam;
    SysSolverCGParams invParam;
    multi1d<int> lattice_size;
    multi1d<LatticeColorMatrix> u;
    multi1d<LatticeColorMatrix> v;

  };

}
#endif
