// $Id: prec_dwf_qprop_array_sse_w.cc,v 1.2 2005-01-06 03:50:48 edwards Exp $
/*! \file
 *  \brief SSE 5D DWF specific quark propagator solver
 */

#include "chromabase.h"
#include "actions/ferm/qprop/prec_dwf_qprop_array_sse_w.h"

#include <sse_dwf_cg.h>

namespace Chroma
{
  namespace SSEDWF
  {
    ///////////////////////////////////////////////////////////////////////////////
    //! Gauge field reader
    double
    gauge_reader(const void *ptr, void *env, 
		 const int latt_coord[Nd], int mu, int row, int col, int reim)
    {
      /* Translate arg */
      multi1d<LatticeColorMatrix>& u = *(multi1d<LatticeColorMatrix>*)ptr;

      // Get node and index
      multi1d<int> coord(Nd);
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
    fermion_reader_rhs(const void *ptr, void *env, 
		       const int latt_coord[5], int color, int spin, int reim)
    {
      /* Translate arg */
      multi1d<LatticeFermion>& psi = *(multi1d<LatticeFermion>*)ptr;
      int Ls1 = psi.size() - 1;

      // Get node and index
      int s = latt_coord[Nd];
      multi1d<int> coord(Nd);
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
	double(psi[Ls1-s].elem(linear).elem(spin).elem(color).real()) : 
	double(psi[Ls1-s].elem(linear).elem(spin).elem(color).imag());

      if (spin >= Ns/2)
	val *= -1;

      return val;
    }


    //! Fermion field reader
    double
    fermion_reader_guess(const void *ptr, void *env, 
			 const int latt_coord[5], int color, int spin, int reim)
    {
      /* Translate arg */
      multi1d<LatticeFermion>& psi = *(multi1d<LatticeFermion>*)ptr;
      int Ls1 = psi.size() - 1;

      // Get node and index
      int s = latt_coord[Nd];
      multi1d<int> coord(Nd);
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
	double(psi[Ls1-s].elem(linear).elem(spin).elem(color).real()) : 
	double(psi[Ls1-s].elem(linear).elem(spin).elem(color).imag());

      if (spin >= Ns/2)
	val *= -1;

      val *= -0.5;

      return val;
    }


    //! Fermion field writer
    void
    fermion_writer_solver(void *ptr, void *env, 
			  const int latt_coord[5], int color, int spin, int reim,
			  double val)
    {
      /* Translate arg */
      multi1d<LatticeFermion>& psi = *(multi1d<LatticeFermion>*)ptr;
      int Ls1 = psi.size() - 1;

      // Get node and index
      int s = latt_coord[Nd];
      multi1d<int> coord(Nd);
      coord = latt_coord;
      int node = Layout::nodeNumber(coord);
      int linear = Layout::linearSiteIndex(coord);

      if (node != Layout::nodeNumber())
      {
	QDPIO::cerr << __func__ << ": wrong coordinates for this node" << endl;
	QDP_abort(1);
      }
 
      // Rescale
      if (spin >= Ns/2)
	val *= -1;

      val *= -2.0;

      // Set the value
      // NOTE: it would be nice to use the "peek" functions, but they will
      // broadcast to all nodes the value since they are platform independent.
      // We don't want that, so we poke into the on-node data
      if (reim == 0)
	psi[Ls1-s].elem(linear).elem(spin).elem(color).real() = val;
      else
	psi[Ls1-s].elem(linear).elem(spin).elem(color).imag() = val;

      return;
    }


    //! Fermion field writer
    void
    fermion_writer_operator(void *ptr, void *env, 
			    const int latt_coord[5], int color, int spin, int reim,
			    double val)
    {
      /* Translate arg */
      multi1d<LatticeFermion>& psi = *(multi1d<LatticeFermion>*)ptr;
      int Ls1 = psi.size() - 1;

      // Get node and index
      int s = latt_coord[Nd];
      multi1d<int> coord(Nd);
      coord = latt_coord;
      int node = Layout::nodeNumber(coord);
      int linear = Layout::linearSiteIndex(coord);

      if (node != Layout::nodeNumber())
      {
	QDPIO::cerr << __func__ << ": wrong coordinates for this node" << endl;
	QDP_abort(1);
      }
 
      // Rescale
      if (spin >= Ns/2)
	val *= -1;

      val *= -0.5;

      // Set the value
      // NOTE: it would be nice to use the "peek" functions, but they will
      // broadcast to all nodes the value since they are platform independent.
      // We don't want that, so we poke into the on-node data
      if (reim == 0)
	psi[Ls1-s].elem(linear).elem(spin).elem(color).real() = val;
      else
	psi[Ls1-s].elem(linear).elem(spin).elem(color).imag() = val;

      return;
    }



    ///////////////////////////////////////////////////////////////////////////////
    void
    solve_cg5(multi1d<LatticeFermion> &solution,    // output
	      const multi1d<LatticeColorMatrix> &U, // input
	      double M5,                            // input
	      double m_f,                           // input
	      const multi1d<LatticeFermion> &rhs,   // input
	      const multi1d<LatticeFermion> &x0,    // input
	      double rsd,                           // input
	      int max_iter,                         // input
	      int &out_iter )                       // output
    {
      // Initialize internal structure of the solver
      //    if (SSE_DWF_init(lattice_size, SSE_DWF_FLOAT, NULL, NULL)) {
      //      error("SSE DWF init() failed");
      //    }

      // Construct shifted gauge field
      multi1d<LatticeColorMatrix> V(Nd);
      for (int i = 0; i < Nd; i++)
	V[i] = shift(U[i], -1, i); // as viewed from the destination

      SSE_DWF_Gauge *g = SSE_DWF_load_gauge(&U, &V, NULL, &SSEDWF::gauge_reader);
      SSE_DWF_Fermion *eta = SSE_DWF_load_fermion(&rhs, NULL, &SSEDWF::fermion_reader_rhs);
      SSE_DWF_Fermion *X0 = SSE_DWF_load_fermion(&x0, NULL, &SSEDWF::fermion_reader_guess);
      SSE_DWF_Fermion *res = SSE_DWF_allocate_fermion();

      QDPIO::cout << "Entering SSE DWF solver: rsd = " << rsd
		  << ", max_iterations = " << max_iter
		  << endl;

      double M_0 = -2*(5.0-M5);
      double out_eps;
      out_eps = 0.0;
      out_iter = 0;
      int status = SSE_DWF_cg_solver(res, &out_eps, &out_iter,
				     g, M_0, m_f, X0, eta, 
				     rsd, max_iter);

      QDPIO::cout << "SSE DWF solver: status = " << status
		  << ", iterations = " << out_iter
		  << ", resulting epsilon = " << out_eps
		  << endl;

      SSE_DWF_save_fermion(&solution, NULL, &SSEDWF::fermion_writer_solver, res);

      SSE_DWF_delete_fermion(res);
      SSE_DWF_delete_fermion(X0);
      SSE_DWF_delete_fermion(eta);
      SSE_DWF_delete_gauge(g);

      // SSE_DWF_fini();
    }

  }  // end namespace SSEDWF


  //! Need a real destructor
  SSEDWFQprop::~SSEDWFQprop() {SSE_DWF_fini();}

  //! Private internal initializer
  void 
  SSEDWFQprop::init()
  {
    QDPIO::cout << "entering SSEEvenOddPrecDWFermActArray::init" << endl;

    if (Nd != 4 || Nc != 3)
    {
      QDPIO::cerr << "SSEEvenOddPrecDWFermActArray: only supports Nd=4 and Nc=3" << endl;
      QDP_abort(1);
    }

    if (invParam.invType != CG_INVERTER)
    {
      QDPIO::cerr << "SSE EvenOddPrecDWF qpropT only supports CG" << endl;
      QDP_abort(1);
    }

    multi1d<int> lattice_size(Nd+1);
    lattice_size[Nd] = N5;
    for(int i=0; i < Nd; ++i)
      lattice_size[i] = Layout::lattSize()[i];

    if (SSE_DWF_init(lattice_size.slice(), SSE_DWF_FLOAT, NULL, NULL) != 0)
    {
      QDPIO::cerr << __func__ << ": error in SSE_DWF_init" << endl;
      QDP_abort(1);
    }

    QDPIO::cout << "exiting SSEEvenOddPrecDWFermActArray::init" << endl;
  }



  //! Solver the linear system
  /*!
   * \param psi      quark propagator ( Modify )
   * \param chi      source ( Read )
   * \return number of CG iterations
   */
  int SSEDWFQprop::operator() (multi1d<LatticeFermion>& psi, const multi1d<LatticeFermion>& chi) const
  {
    QDPIO::cout << "entering SSEEvenOddPrecDWFermActArray::qpropT" << endl;

    START_CODE();

    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    // Apply SSE inverter
    double M5  = toDouble(OverMass);
    double m_f = toDouble(Mass);
    double rsd = toDouble(invParam.RsdCG);
    double rsd_sq = rsd * rsd;
    int    max_iter = invParam.MaxCG;
    int    n_count;
    SSEDWF::solve_cg5(psi, u, M5, m_f, chi, psi, rsd_sq, max_iter, n_count);

    END_CODE();

    QDPIO::cout << "exiting SSEEvenOddPrecDWFermActArray::qpropT" << endl;

    return n_count;
  }

}
