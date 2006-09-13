// $Id: prec_dwf_qprop_array_cg_dwf_w.cc,v 3.1 2006-09-13 19:06:17 bjoo Exp $
/*! \file
 *  \brief SSE 5D DWF specific quark propagator solver
 */

#include "chromabase.h"
#include "actions/ferm/qprop/prec_dwf_qprop_array_cg_dwf_w.h"

/* Always use the double prec solver */
/* Could reuse almost exactly the same code for altivec and bgl
   but annoyingly the silly names are defined by the L3 macro. I'll see
   if I can work around this later. For now do it all in DP */
#include <dwf-ssed.h>

namespace Chroma
{
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
				    g, M_0, m_f, X0, eta, 
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


  //----------------------------------------------------------------------------------
  //! Private internal initializer
  void CGDWFQpropT::init(Handle< FermState<LatticeFermion,
			  multi1d<LatticeColorMatrix>,
			  multi1d<LatticeColorMatrix> > > state,
			  const GroupXML_t& inv)
  {
    QDPIO::cout << "entering CGDWFQpropT::init" << endl;

    if (Nd != 4 || Nc != 3)
    {
      QDPIO::cerr << "CGDWFQpropT: only supports Nd=4 and Nc=3" << endl;
      QDP_abort(1);
    }

    // Read the XML for the CG params
    try
    {
      std::istringstream  is(inv.xml);
      XMLReader  paramtop(is);

      read(paramtop, inv.path, invParam);
    }
    catch (const std::string& e)
    {
      QDPIO::cerr << "CGDWFQpropT: only support a CG inverter" << endl;
      QDP_abort(1);
    }

    multi1d<int> lattice_size(Nd+1);
    lattice_size[Nd] = N5;
    for(int i=0; i < Nd; ++i)
      lattice_size[i] = Layout::lattSize()[i];

    if (N5 % 4 != 0)
    {
      QDPIO::cerr << "SSE qpropT only N5 that is multiple of 2" << endl;
      QDP_abort(1);
    }

    std::string dwf_error_str = "double prec.";

    if (MIT_ssed_DWF_init(lattice_size.slice(), NULL, NULL) != 0)
    {
      QDPIO::cerr << __func__ << ": error in MIT_ssed_DWF_init: " << dwf_error_str << endl;
      QDP_abort(1);
    }

    // Transform the U fields
    multi1d<LatticeColorMatrix> u = state->getLinks();
  
    if (anisoParam.anisoP)
    {
      Real ff = where(anisoParam.anisoP, anisoParam.nu / anisoParam.xi_0, Real(1));

      // Rescale the u fields by the anisotropy
      for(int mu=0; mu < u.size(); ++mu)
      {
	if (mu != anisoParam.t_dir)
	  u[mu] *= ff;
      }
    }

    // Construct shifted gauge field
    multi1d<LatticeColorMatrix> v(Nd);
    for (int i = 0; i < Nd; i++)
      v[i] = shift(u[i], -1, i); // as viewed from the destination

    g = MIT_ssed_DWF_load_gauge(&u, &v, NULL, &CGDWF::gauge_reader);

    QDPIO::cout << "exiting CGDWFQpropT::init" << endl;
  }


  //----------------------------------------------------------------------------------
  //! Need a real destructor
  void CGDWFQpropT::fini()
  {
    QDPIO::cout << "CGDWFQpropT: calling destructor" << endl;
    MIT_ssed_DWF_delete_gauge(g);
    MIT_ssed_DWF_fini();
  }


  //----------------------------------------------------------------------------------
  //! Solver the linear system
  /*!
   * \param psi      quark propagator ( Modify )
   * \param chi      source ( Read )
   * \return number of CG iterations
   */
  SystemSolverResults_t 
  CGDWFQpropT::operator() (multi1d<LatticeFermion>& psi, const multi1d<LatticeFermion>& chi) const
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
    CGDWF::solve_cg5(psi, g, M5, m_f, 
		      chi, psi, rsd_sq, max_iter, out_eps, res.n_count);

//    fini();   // only needed because 2 qpropT might be active - SSE CG does not allow this

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

}
