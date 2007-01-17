// -*- C++ -*-
// $Id: avp_inverter_interface.h,v 3.8 2007-01-17 02:39:27 edwards Exp $
/*! \file
 *  \brief Base class for AVP's DWF solver interface
 */

#ifndef AVP_INVERTER_INTERFACE_H
#define AVP_INVERTER_INTERFACE_H

#include "chromabase.h"

using namespace QDP;

namespace Chroma 
{ 
  //! AVP's DWF Solver interface
  /*!
   * \ingroup qprop
   *
   * @{
   */
  namespace AVPSolverFunctions 
  { 
    // Gauge Reader function - user supplied
    double gaugeReader(const void *OuterGauge,
		       void *env,
		       const int lattice_coord[4],
		       int mu,
		       int row,
		       int col,
		       int reim);
      
      
    // Fermion Reader function - user supplied
    double fermionReaderRHS(const void *OuterFermion,
			    void *env,
			    const int lattice_coord[5],
			    int color,
			    int spin,
			    int reim);
      
    // Fermion Reader function - user supplied
    double fermionReaderGuess(const void *OuterFermion,
			      void *env,
			      const int lattice_addr[5],
			      int color,
			      int spin,
			      int reim) ;
      
      
    // Fermion Writer function - user supplied
    void fermionWriterSolver(void *OuterFermion, 
			     void *env, 
			     const int lattice_addr[5],
			     int color, 
			     int spin,
			     int reim,
			     double val);
	
      
    // Fermion Writer function - user supplied
    void fermionWriterOperator(void *OuterFermion, 
			       void *env, 
			       const int lattice_addr[5],
			       int color, 
			       int spin,
			       int reim,
			       double val);

  }  // namespace AVPSolverFunctions 

  namespace AVPSolver 
  {

    template< typename U, typename T>
    class AVPSolverInterface 
    {
    public:
      
      // Gauge Field Type
      typedef U Gauge;
      
      // Gauge Fermion Field Type
      typedef T Fermion;
      
      
    protected:
      
      // Here I can wrap AVPs solver stuff
      virtual Fermion* loadFermionRHS(const void *OuterFermion) const = 0;
      virtual Fermion* loadFermionGuess(const void *OuterFermion) const = 0; 
      virtual Fermion* allocateFermion(void) const =0 ;
      virtual void saveFermionSolver(void *OuterFermion,
				     Fermion* CGfermion) const =0;

      virtual void saveFermionOperator(void *OuterFermion,
				       Fermion* CGfermion) const =0;
      
      virtual void deleteFermion(Fermion* ptr) const = 0;
      
      
      
      virtual int cgInternal(Fermion       *psi,
			     double        *out_eps,
			     int           *out_iter,
			     double        M,
			     double        m_f,
			     const Fermion *x0,
			     const Fermion *eta,
			     double        eps,
			     int           min_iter,
			     int           max_iter) const = 0;



      

      
    public:

      // Init the system -- Constructor call?
      virtual int init(const int lattice[5],
		       void *(*allocator)(size_t size),
		       void (*deallocator)(void *)) = 0;
      
      // Finalize - destructor call
      virtual void fini(void) = 0;
      
      virtual void loadGauge(const void *OuterGauge_U,
			     const void *OuterGauge_V) = 0;
      
      virtual void deleteGauge(void) = 0;


      // Call the solver
      int cgSolver(multi1d<LatticeFermion> &solution,    // output
		   double M5,                            // input
		   double m_f,                           // input
		   const multi1d<LatticeFermion> &rhs,   // input
		   const multi1d<LatticeFermion> &x0,    // input
		   double rsd,                           // input
		   int max_iter,                         // input
		   double& out_eps,                      // output
		   int &out_iter )       const           // output
	{

	  StopWatch swatch;
	  swatch.reset();
	  swatch.start();
	  Fermion *eta = loadFermionRHS(&rhs);
	  swatch.stop();
	  QDPIO::cout << "Importing RHS Fermion took: " << swatch.getTimeInSeconds() << " seconds " << endl;


	  swatch.reset();
	  swatch.start();
	  Fermion *X0  = loadFermionGuess(&x0);
	  swatch.stop();

	  QDPIO::cout << "Importing Guess took: " << swatch.getTimeInSeconds() 
		      << " seconds" << endl;
	
	  Fermion *res = allocateFermion();
	
	  QDPIO::cout << "Entering CG_DWF solver: rsd = " << rsd
		      << ", max_iterations = " << max_iter
		      << endl;
	
	  double M_0 = -2*M5;
	  out_eps = 0.0;
	  out_iter = 0;
	  int min_iter = 0;


	  swatch.reset();
	  swatch.start();
	
	
	  int status = cgInternal(res, &out_eps, &out_iter,
				  M_0, m_f, X0, eta, 
				  rsd, min_iter, max_iter);
	
	  swatch.stop();
	  QDPIO::cout << "CGInternal : status = " << status
		      << ", iterations = " << out_iter
		      << ", resulting epsilon = " << out_eps
		      << endl;
	
	  if (status != 0) {
	  
	    QDPIO::cerr << "DWF_solver: status = " << status
			<< ", iterations = " << out_iter
			<< ", resulting epsilon = " << out_eps
			<< endl;
	    QDP_abort(1);
	  }

	  // Flop counting
	  {
	    unsigned long Ls = solution.size();
	    unsigned long Ndiag  = (4*Ls+2)*Nc*Ns; /* This is my count with the blas / chiral proj ops */
	    unsigned long NdiagInv = (10*Ls-8)*Nc*Ns;
	    unsigned long Neo    = Ls*(1320+24);
	    unsigned long N_mpsi = 2*Ndiag + 2*Neo + Ls*24;
	    unsigned long Nflops_c = (24*Ls + 2*N_mpsi) + (48*Ls);  /* constant term */
	    unsigned long Nflops_s = (2*N_mpsi + Ls*(2*48+2*24));   /* slope */
	    unsigned long long Nflops_per_cbsite = Nflops_c + ( out_iter)*Nflops_s;
	    unsigned long long Nflops_total = Nflops_per_cbsite*(Layout::sitesOnNode()/2);
	  
	    /* Flop count for inverter */
	    FlopCounter flopcount;
	    flopcount.reset();
	    flopcount.addFlops(Nflops_total);
	    flopcount.report("CGDWFQpropT", swatch.getTimeInSeconds());
	  }

	  swatch.reset();
	  swatch.start();
	  saveFermionSolver(&solution,res);
	  swatch.stop();

	  QDPIO::cout << "Exporting Solution took: " << swatch.getTimeInSeconds()
		      << " seconds " << endl;

	  deleteFermion(res);
	  deleteFermion(X0);
	  deleteFermion(eta);
	}
      
    }; // Class
  } // AVP Solver Interface

  /*! @} */   // end of group qprop

} // Chroma namepace

#endif
