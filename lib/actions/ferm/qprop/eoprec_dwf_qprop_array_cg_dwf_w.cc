// $Id: eoprec_dwf_qprop_array_cg_dwf_w.cc,v 3.1 2006-10-19 16:01:33 edwards Exp $
/*! \file
 *  \brief SSE 5D DWF specific quark propagator solver
 */

#include "chromabase.h"
#include "actions/ferm/qprop/eoprec_dwf_qprop_array_cg_dwf_w.h"

/* Always use the double prec solver */
/* Could reuse almost exactly the same code for altivec and bgl
   but annoyingly the silly names are defined by the L3 macro. I'll see
   if I can work around this later. For now do it all in DP */
#include <dwf-ssed.h>

namespace Chroma
{

#if 0
  namespace CGDWF
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
	      MIT_ssed_DWF_Gauge* g,                     // input
	      double M5,                            // input
	      double m_f,                           // input
	      const multi1d<LatticeFermion> &rhs,   // input
	      const multi1d<LatticeFermion> &x0,    // input
	      double rsd,                           // input
	      int max_iter,                         // input
	      double& out_eps,                      // output
	      int &out_iter )                       // output
    {

      MIT_ssed_DWF_Fermion *eta = MIT_ssed_DWF_load_fermion(&rhs, NULL, &CGDWF::fermion_reader_rhs);
      MIT_ssed_DWF_Fermion *X0 = MIT_ssed_DWF_load_fermion(&x0, NULL, &CGDWF::fermion_reader_guess);
      MIT_ssed_DWF_Fermion *res = MIT_ssed_DWF_allocate_fermion();

      QDPIO::cout << "Entering MIT_ssed_DWF solver: rsd = " << rsd
		  << ", max_iterations = " << max_iter
		  << endl;

      double M_0 = -2*M5;
      out_eps = 0.0;
      out_iter = 0;
      int min_iter=5;

      StopWatch swatch;
      swatch.reset();
      swatch.start();

      int status = MIT_ssed_DWF_cg_solver(res, &out_eps, &out_iter,
				    (SinglePrecSolver::Gauge *)g, 
				    M_0, m_f, X0, eta, 
				    rsd, min_iter, max_iter);

      swatch.stop();
      QDPIO::cout << "MIT_ssed_DWF_cg_solver: status = " << status
		  << ", iterations = " << out_iter
		  << ", resulting epsilon = " << out_eps
                  << endl;

      if (status != 0)
      {
	QDPIO::cerr << "Error in MIT_ssed_DWF_solver: status = " << status
		    << ", iterations = " << out_iter
		    << ", resulting epsilon = " << out_eps
		    << endl;
	QDP_abort(1);
      }

      // Flop counting
      {
	unsigned long Ls = x0.size();
	unsigned long Ndiag  = (4*Ls+2)*Nc*Ns; /* This is my count with the blas / chiral proj ops */
	unsigned long NdiagInv = (10*Ls-8)*Nc*Ns;
	unsigned long Neo    = Ls*(1320+24);
	unsigned long N_mpsi = 2*Ndiag + 2*Neo + Ls*24;
	unsigned long Nflops_c = (24*Ls + 2*N_mpsi) + (48*Ls);  /* constant term */
	unsigned long Nflops_s = (2*N_mpsi + Ls*(2*48+2*24));   /* slope */
	unsigned long long Nflops_per_cbsite = Nflops_c + out_iter*Nflops_s;
	unsigned long long Nflops_total = Nflops_per_cbsite*(Layout::sitesOnNode()/2);

	/* Flop count for inverter */
	FlopCounter flopcount;
	flopcount.reset();
	flopcount.addFlops(Nflops_total);
	flopcount.report("CGDWFQpropT", swatch.getTimeInSeconds());
      }

      MIT_ssed_DWF_save_fermion(&solution, NULL, &CGDWF::fermion_writer_solver, res);

      MIT_ssed_DWF_delete_fermion(res);
      MIT_ssed_DWF_delete_fermion(X0);
      MIT_ssed_DWF_delete_fermion(eta);

      // SSE_DWF_fini();
    }

  }  // end namespace CGDWF
#endif

  //----------------------------------------------------------------------------------
  //! Private internal initializer
  
}
