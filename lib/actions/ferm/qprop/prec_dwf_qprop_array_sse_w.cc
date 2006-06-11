// $Id: prec_dwf_qprop_array_sse_w.cc,v 3.3 2006-06-11 18:45:35 edwards Exp $
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
	      SSE_DWF_Gauge* g,                     // input
	      double M5,                            // input
	      double m_f,                           // input
	      const multi1d<LatticeFermion> &rhs,   // input
	      const multi1d<LatticeFermion> &x0,    // input
	      double rsd,                           // input
	      int max_iter,                         // input
	      double& out_eps,                      // output
	      int &out_iter )                       // output
    {
      // Initialize internal structure of the solver
      //    if (SSE_DWF_init(lattice_size, SSE_DWF_FLOAT, NULL, NULL)) {
      //      error("SSE DWF init() failed");
      //    }

      SSE_DWF_Fermion *eta = SSE_DWF_load_fermion(&rhs, NULL, &SSEDWF::fermion_reader_rhs);
      SSE_DWF_Fermion *X0 = SSE_DWF_load_fermion(&x0, NULL, &SSEDWF::fermion_reader_guess);
      SSE_DWF_Fermion *res = SSE_DWF_allocate_fermion();

      QDPIO::cout << "Entering SSE DWF solver: rsd = " << rsd
		  << ", max_iterations = " << max_iter
		  << endl;

      double M_0 = -2*M5;
      out_eps = 0.0;
      out_iter = 0;

      StopWatch swatch;
      swatch.reset();
      swatch.start();

      int status = SSE_DWF_cg_solver(res, &out_eps, &out_iter,
				     g, M_0, m_f, X0, eta, 
				     rsd, max_iter);

      swatch.stop();
      QDPIO::cout << "SSE DWF solver: status = " << status
		  << ", iterations = " << out_iter
		  << ", resulting epsilon = " << out_eps
                  << endl;

      if (status != 0)
      {
	QDPIO::cerr << "Error in SSE DWF solver: status = " << status
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
	flopcount.report("SSEDWFQpropT", swatch.getTimeInSeconds());
      }

      SSE_DWF_save_fermion(&solution, NULL, &SSEDWF::fermion_writer_solver, res);

      SSE_DWF_delete_fermion(res);
      SSE_DWF_delete_fermion(X0);
      SSE_DWF_delete_fermion(eta);

      // SSE_DWF_fini();
    }

  }  // end namespace SSEDWF


  //----------------------------------------------------------------------------------
  //! Private internal initializer
  void SSEDWFQpropT::init(Handle< FermState<LatticeFermion,
			  multi1d<LatticeColorMatrix>,
			  multi1d<LatticeColorMatrix> > > state)
  {
    QDPIO::cout << "entering SSEDWFQpropT::init" << endl;

    if (Nd != 4 || Nc != 3)
    {
      QDPIO::cerr << "SSEDWFQpropT: only supports Nd=4 and Nc=3" << endl;
      QDP_abort(1);
    }

    if (invParam.invType != CG_INVERTER)
    {
      QDPIO::cerr << "SSE qpropT only supports CG" << endl;
      QDP_abort(1);
    }

    multi1d<int> lattice_size(Nd+1);
    lattice_size[Nd] = N5;
    for(int i=0; i < Nd; ++i)
      lattice_size[i] = Layout::lattSize()[i];

    if (N5 % 4 != 0)
    {
      QDPIO::cerr << "SSE qpropT only N5 that is multiple of 4" << endl;
      QDP_abort(1);
    }

#if BASE_PRECISION == 32
    SSE_DWF_FP_SIZE dwf_prec = SSE_DWF_FLOAT;
    std::string dwf_error_str = "single prec.";
#elif BASE_PRECISION == 64
    SSE_DWF_FP_SIZE dwf_prec = SSE_DWF_DOUBLE;
    std::string dwf_error_str = "double prec.";
#else
#error "Unknown BASE_PRECISION"
#endif

    if (SSE_DWF_init(lattice_size.slice(), dwf_prec, NULL, NULL) != 0)
    {
      QDPIO::cerr << __func__ << ": error in SSE_DWF_init: " << dwf_error_str << endl;
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

    g = SSE_DWF_load_gauge(&u, &v, NULL, &SSEDWF::gauge_reader);

    QDPIO::cout << "exiting SSEDWFQpropT::init" << endl;
  }


  //----------------------------------------------------------------------------------
  //! Need a real destructor
  void SSEDWFQpropT::fini()
  {
    QDPIO::cout << "SSEDWFQpropT: calling destructor" << endl;
    SSE_DWF_delete_gauge(g);
    SSE_DWF_fini();
  }


  //----------------------------------------------------------------------------------
  //! Solver the linear system
  /*!
   * \param psi      quark propagator ( Modify )
   * \param chi      source ( Read )
   * \return number of CG iterations
   */
  SystemSolverResults_t 
  SSEDWFQpropT::operator() (multi1d<LatticeFermion>& psi, const multi1d<LatticeFermion>& chi) const
  {
    QDPIO::cout << "entering SSEDWFQpropT::operator()" << endl;

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
    SSEDWF::solve_cg5(psi, g, M5, m_f, 
		      chi, psi, rsd_sq, max_iter, out_eps, res.n_count);

//    fini();   // only needed because 2 qpropT might be active - SSE CG does not allow this

    res.resid = out_eps;

    QDPIO::cout << "exiting SSEDWFQpropT::operator()" << endl;

    END_CODE();

    return res;
  }

}
