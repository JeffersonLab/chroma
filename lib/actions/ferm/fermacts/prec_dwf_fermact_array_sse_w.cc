// $Id: prec_dwf_fermact_array_sse_w.cc,v 1.2 2004-09-06 16:16:24 edwards Exp $
/*! \file
 *  \brief SSE 4D style even-odd preconditioned domain-wall fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_array_sse_w.h"
#include "actions/ferm/linop/unprec_dwf_linop_array_w.h"
#include "actions/ferm/linop/prec_dwf_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include <sse_dwf_cg.h>

//! Private internal initializer
void 
SSEEvenOddPrecDWFermActArray::init()
{
  QDPIO::cout << "entering SSEEvenOddPrecDWFermActArray::init" << endl;

  if (Nd != 4 || Nc != 3)
  {
    QDPIO::cerr << "SSEEvenOddPrecDWFermActArray: only supports Nd=4 and Nc=3" << endl;
    QDP_abort(1);
  }

  multi1d<int> lattice_size(Nd+1);
  lattice_size[Nd] = size();
  for(int i=0; i < Nd; ++i)
    lattice_size[i] = Layout::lattSize()[i];

  if (SSE_DWF_init(lattice_size.slice(), SSE_DWF_FLOAT, NULL, NULL) != 0)
  {
    QDPIO::cerr << __func__ << ": error in SSE_DWF_init" << endl;
    QDP_abort(1);
  }

  QDPIO::cout << "exiting SSEEvenOddPrecDWFermActArray::init" << endl;
}


//! Private internal destructor
void 
SSEEvenOddPrecDWFermActArray::fini()
{
  SSE_DWF_fini();
}



//! Produce a linear operator for this action
/*!
 * \ingroup fermact
 *
 * The operator acts on the odd sublattice
 *
 * \param state 	    gauge field     	       (Read)
 */
const EvenOddPrecLinearOperator<multi1d<LatticeFermion> >*
SSEEvenOddPrecDWFermActArray::linOp(Handle<const ConnectState> state) const
{
  return new EvenOddPrecDWLinOpArray(state->getLinks(),WilsonMass,m_q,N5);
}

//! Produce a M^dag.M linear operator for this action
/*!
 * The operator acts on the odd sublattice
 *
 * \param state 	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >*
SSEEvenOddPrecDWFermActArray::lMdagM(Handle<const ConnectState> state) const
{
  return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
}

//! Produce a linear operator for this action but with quark mass 1
/*!
 * The operator acts on the entire lattice
 *
 * \param state	    gauge field     	       (Read)
 */
const LinearOperator<multi1d<LatticeFermion> >*
SSEEvenOddPrecDWFermActArray::linOpPV(Handle<const ConnectState> state) const
{
  // For the PV operator, use the **unpreconditioned** one
  // fixed to quark mass 1
  return new UnprecDWLinOpArray(state->getLinks(),WilsonMass,1.0,N5);
}


//! Gauge field reader
double
DWF_gauge_reader(const void *ptr, void *env, 
		 const int latt_coord[4], int mu, int row, int col, int reim)
{
  /* Translate arg */
  multi1d<LatticeColorMatrix>& u = *(multi1d<LatticeColorMatrix>*)ptr;

  // Get node and index
  multi1d<int> coord(4);
  coord = latt_coord;
  int node = Layout::nodeNumber(coord);
  int linear = Layout::linearSiteIndex(coord);

  if (node != Layout::nodeNumber())
  {
    QDPIO::cerr << __func__ << ": wrong coordinates for this node" << endl;
    QDP_abort(1);
  }
 
  // Get the value
  // NOTE: it would be nice to use the "peek" functions, but they will
  // broadcast to all nodes the value since they are platform independent.
  // We don't want that, so we poke into the on-node data
  double val = (reim == 0) ? 
    toDouble(u[mu].elem(linear).elem().elem(row,col).real()) : 
    toDouble(u[mu].elem(linear).elem().elem(row,col).imag());

  return val;
}


//! Fermion field reader
double
DWF_fermion_reader(const void *ptr, void *env, 
		   const int latt_coord[5], int color, int spin, int reim)
{
  /* Translate arg */
  multi1d<LatticeFermion>& psi = *(multi1d<LatticeFermion>*)ptr;

  // Get node and index
  int s = latt_coord[4];
  multi1d<int> coord(4);
  coord = latt_coord;
  int node = Layout::nodeNumber(coord);
  int linear = Layout::linearSiteIndex(coord);

  if (node != Layout::nodeNumber())
  {
    QDPIO::cerr << __func__ << ": wrong coordinates for this node" << endl;
    QDP_abort(1);
  }
 
  // Get the value
  // NOTE: it would be nice to use the "peek" functions, but they will
  // broadcast to all nodes the value since they are platform independent.
  // We don't want that, so we poke into the on-node data
  double val = (reim == 0) ? 
    double(psi[s].elem(linear).elem(spin).elem(color).real()) : 
    double(psi[s].elem(linear).elem(spin).elem(color).imag());

  return val;
}


//! Fermion field reader
double
DWF_fermion_reader_zero(const void *ptr, void *env, 
			const int latt_coord[5], int color, int spin, int reim)
{
  return 0.0;
}


//! Fermion field reader
void
DWF_fermion_writer(void *ptr, void *env, 
		   const int latt_coord[5], int color, int spin, int reim,
		   double val)
{
  /* Translate arg */
  multi1d<LatticeFermion>& psi = *(multi1d<LatticeFermion>*)ptr;

  // Get node and index
  int s = latt_coord[4];
  multi1d<int> coord(4);
  coord = latt_coord;
  int node = Layout::nodeNumber(coord);
  int linear = Layout::linearSiteIndex(coord);

  if (node != Layout::nodeNumber())
  {
    QDPIO::cerr << __func__ << ": wrong coordinates for this node" << endl;
    QDP_abort(1);
  }
 
  // Set the value
  // NOTE: it would be nice to use the "peek" functions, but they will
  // broadcast to all nodes the value since they are platform independent.
  // We don't want that, so we poke into the on-node data
  if (reim == 0)
    psi[s].elem(linear).elem(spin).elem(color).real() = val;
  else
    psi[s].elem(linear).elem(spin).elem(color).imag() = val;

  return;
}



//! Optimized inverter - this is temporary
void 
SSEEvenOddPrecDWFermActArray::opt_qpropT(multi1d<LatticeFermion>& psi, 
					 Handle<const ConnectState> state, 
					 const multi1d<LatticeFermion>& chi, 
					 enum InvType invType,
					 const Real& RsdCG, 
					 int MaxCG, int& ncg_had) const
{
  QDPIO::cout << "entering SSEEvenOddPrecDWFermActArray::opt_qpropT" << endl;

  START_CODE();

  const multi1d<LatticeColorMatrix>& u = state->getLinks();

  // Get a shift copy of the gauge fields
  multi1d<LatticeColorMatrix> v(u.size());
  for(int mu=0; mu < Nd; ++mu)
    v[mu] = shift(u[mu], BACKWARD, mu);

  // Initialize the SSE specific fields
  SSE_DWF_Gauge*   g   = SSE_DWF_load_gauge(&u, &v, NULL, DWF_gauge_reader);
  if (g == NULL)
    QMP_error_exit("%s: error allocating g", __func__);

  multi1d<LatticeFermion> chi_tmp(N5);
  for(int s=0; s < N5; ++s)
  {
//    chi_tmp[s] = -Real(2)*chi[s];   // compensate for AVP's def
    chi_tmp[s] = chi[s];
  }
 
  SSE_DWF_Fermion* rhs = SSE_DWF_load_fermion(&chi_tmp, NULL, DWF_fermion_reader);
  if (rhs == NULL)
    QMP_error_exit("%s: error allocating rhs", __func__);

  SSE_DWF_Fermion* x0  = SSE_DWF_load_fermion(&psi, NULL, DWF_fermion_reader);
  if (x0 == NULL)
    QMP_error_exit("%s: error allocating x0", __func__);

  SSE_DWF_Fermion* x   = SSE_DWF_allocate_fermion();
  if (x == NULL)
    QMP_error_exit("%s: error allocating x", __func__);


  // Call the solver
  double out_epsilon = 0.0;
  double M_0 = -2*toDouble(Double(a5*(Nd-WilsonMass) + 1));  // compensate for AVP's def
  double mq  = toDouble(m_q);
  double rsdsq = toDouble(RsdCG)*toDouble(RsdCG);
  int err;

  cerr << "M_0=" << M_0 << " mq=" << mq << " rsdsq=" << rsdsq << endl;

 
  err = SSE_DWF_cg_solver(x, &out_epsilon, &ncg_had,
			  g, M_0, mq,
			  x0, rhs, 
			  rsdsq, MaxCG);

  QDPIO::cerr << __func__ 
	      << ": iterations = " << ncg_had 
	      << "  eps = " << sqrt(out_epsilon)
	      << endl;

  if (err != 0)
  {
    QDPIO::cerr << __func__ << ": convergence not found" << endl;
    QDP_abort(1);
  }

  // Save the result
  psi = zero;
  SSE_DWF_save_fermion(&psi, NULL, DWF_fermion_writer, x);

  // Cleanup
  SSE_DWF_delete_fermion(x);
  SSE_DWF_delete_fermion(x0);
  SSE_DWF_delete_fermion(rhs);
  SSE_DWF_delete_gauge(g);

  END_CODE();

  QDPIO::cout << "exiting SSEEvenOddPrecDWFermActArray::opt_qpropT" << endl;
}
