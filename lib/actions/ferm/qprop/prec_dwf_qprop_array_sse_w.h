// -*- C++ -*-
// $Id: prec_dwf_qprop_array_sse_w.h,v 1.1 2005-01-02 05:21:10 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#ifndef __prec_dwf_qprop_array_sse_w_h__
#define __prec_dwf_qprop_array_sse_w_h__

#include "fermact.h"

using namespace QDP;

namespace Chroma
{

  //! SSE Propagator DWF qpropT
  /*! \ingroup qprop
   *
   * Propagator solver for DWF fermions
   */
  class SSEDWFQprop : public SystemSolver< multi1d<LatticeFermion> >
  {
  public:
    //! Constructor
    /*!
     * \param qpropT_    5D solver ( Read )
     * \param PV_        Pauli-Villars linear operator ( Read )
     * \param m_q_       quark mass ( Read )
     */
    SSEDWFQprop(Handle<const ConnectState> state_, 
		const Real& OverMass,
		const Real& Mass,
		const int N5,
		const InvertParam_t& invParam_) : 
      state(state_), OverMass(OverMass_), Mass(Mass_), N5(N5_), invParam(invParam_) {init();}

    //! Need a real destructor
    ~SSEDWFQprop();

    //! Return the subset on which the operator acts
    const OrderedSubset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    int operator() (LatticeFermion& psi, const LatticeFermion& chi) const;

  protected:
    //! Private internal initializer
    void init();
    double gauge_reader(const void *ptr, void *env, 
			const int latt_coord[Nd], int mu, int row, int col, int reim);

    double fermion_reader_rhs(const void *ptr, void *env, 
			      const int latt_coord[5], int color, int spin, int reim);

    double fermion_reader_guess(const void *ptr, void *env, 
				const int latt_coord[5], int color, int spin, int reim);

    void fermion_writer_solver(void *ptr, void *env, 
			       const int latt_coord[5], int color, int spin, int reim,
			       double val);

    void fermion_writer_operator(void *ptr, void *env, 
				 const int latt_coord[5], int color, int spin, int reim,
				 double val);

    void solve_cg5(multi1d<LatticeFermion> &solution,    // output
		   const multi1d<LatticeColorMatrix> &U, // input
		   double M5,                            // input
		   double m_f,                           // input
		   const multi1d<LatticeFermion> &rhs,   // input
		   const multi1d<LatticeFermion> &x0,    // input
		   double rsd,                           // input
		   int max_iter,                         // input
		   int &out_iter );                      // output
      
  private:
    // Hide default constructor
    SSEDWFQprop() {}

    Handle<const ConnectState> state;
    const Real OverMass;
    const Real Mass;
    const int  N5;
    const InvertParam_t invParam;
  };

}

#endif
