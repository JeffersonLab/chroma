// -*- C++ -*-
// $Id: mdwf_solver.h,v 1.1 2008-05-05 19:18:06 bjoo Exp $
/*! \file
 *  \brief DWF/SSE double-prec solver
 */

#ifndef MDWF_SOLVER_H
#define MDWF_SOLVER_H


extern "C" {
#include <qop-mdwf3.h>
};

#include "eoprec_constdet_wilstype_fermact_w.h"
#include "io/aniso_io.h"
#include "actions/ferm/invert/syssolver_cg_params.h"


using namespace QDP;
namespace Chroma 
{ 
  //! AVP's DWF Solver interface
  /*!
   * \ingroup qprop
   *
   * @{
   */
  namespace { 
    //      Gauge Reader function - user supplied
    double gaugeReader(int mu,	
		       const int latt_coord[4],
		       int row,
		       int col,
		       int reim,
		       void *env)
    {
	/* Translate arg */
      multi1d<LatticeColorMatrix>& u = *(multi1d<LatticeColorMatrix>*)env;
      
      // Get node and index
      multi1d<int> coord(Nd);
      coord = latt_coord;
      int node = Layout::nodeNumber(coord);
      int linear = Layout::linearSiteIndex(coord);
      
      if (node != Layout::nodeNumber()) {
	
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
    
    
    // Fermion Reader function - user supplied
    double fermionReader(const int latt_coord[5],
			 int color,
			 int spin,
			 int reim,
			 void *env)
    {
      /* Translate arg */
      multi1d<LatticeFermion>& psi = *(multi1d<LatticeFermion>*)env;
      
      // Get node and index
      int s = latt_coord[Nd];
      multi1d<int> coord(Nd);
      coord = latt_coord;
      int node = Layout::nodeNumber(coord);
      int linear = Layout::linearSiteIndex(coord);
      
      if (node != Layout::nodeNumber()) {
	
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
      
      // if (spin >= Ns/2)
      //	val *= -1;
      
      return val;
    }
    
      
      
    
    // Fermion Writer function - user supplied
    void fermionWriter(const int latt_coord[5],
		       int color, 
		       int spin,
		       int reim,
		       double val,
		       void *env )
      
    {
      /* Translate arg */
      multi1d<LatticeFermion>& psi = *(multi1d<LatticeFermion>*)env;
      
      // Get node and index
      int s = latt_coord[Nd];
      multi1d<int> coord(Nd);
      coord = latt_coord;
      int node = Layout::nodeNumber(coord);
      int linear = Layout::linearSiteIndex(coord);
      
      if (node != Layout::nodeNumber()) {
	QDPIO::cerr << __func__ << ": wrong coordinates for this node" << endl;
	QDP_abort(1);
      }
      
      // Rescale
      //      if (spin >= Ns/2)
      //	val *= -1;
      
      // val *= -2.0;
      
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
    
    // Env is for the interface spec. I ignore it completely
    void sublattice_func(int lo[],
			 int hi[],
			 const int node[],
			 void *env) {
      
      // Given my node coordinates in node[]
      // produce the lo/hi pair
      
      // This is the size on the local subgrid. 
      // For QDP++ they are all the same.
      const multi1d<int>& local_subgrid=Layout::subgridLattSize();
      for(int i=0; i <Nd; i++) { 
	// The lowest coordinate is just the node coordinate
	// times the local subgrid size in that direction
	lo[i]=node[i]*local_subgrid[i];
	
	// The high is the start of the next 'corner'
	// I would say my hi is hi[i]-1 but am following the 
	// document conventions
	hi[i]=(node[i]+1)*local_subgrid[i];
      }
    }
  }

  class MDWFQpropT : public SystemSolverArray<LatticeFermion>  { 
    
  public:
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;
    
    /* Constructor */
    MDWFQpropT(Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A_,
	       Handle< LinOpSystemSolverArray<T> > invA_,   // throw away
	       Handle< FermState<T,P,Q> > state_, 
	       const Real& OverMass_,
	       const Real& Mass_,
	       const AnisoParam_t& anisoParam_,
	       const GroupXML_t& invParam_) : A(A_), 
					      OverMass(OverMass_), 
					      Mass(Mass_),       
					      N5(A->size()), 
					      anisoParam(anisoParam_)  {
           init(state_, invParam_);
    } 
    
      
    /* Destructor */  
    ~MDWFQpropT() { 
      fini(); state = NULL; 
    }
      
    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    SystemSolverResults_t operator() (multi1d<LatticeFermion>& psi, 
				      const multi1d<LatticeFermion>& chi) const
    {
      SystemSolverResults_t res;
      int out_iters_single=0;
      int out_iters_double=0;

      /* Single Precision Branch. Test this first
	 then duplicate for double prec */
#if 0
      {
	double out_eps_single;
	
	// OK Try to solve for the single precision 
	// Cast to single precision gauge field
	QOP_F3_MDWF_Gauge *sprec_gauge = NULL;
	if ( QOP_F3_MDWF_import_gauge(&sprec_gauge,
				      state,
				      gaugeReader,
				      (void *)&u) != 0 ) { 
	  QDPIO::cerr << "MDWF Error:  "<< QOP_MDWF_error(state) << endl;
	  QDP_abort(1);
	}
	
	/* Fermion fields */
	QOP_F3_MDWF_Fermion *sprec_rhs;
	QOP_F3_MDWF_Fermion *sprec_x0;
	QOP_F3_MDWF_Fermion *sprec_soln;
	
	if( QOP_F3_MDWF_import_fermion(&sprec_rhs, 
				       state,
				       fermionReader,
				       (void *)&chi) != 0 ) {
	  QDPIO::cerr << "MDWF Error:  "<< QOP_MDWF_error(state) << endl;
	  QDP_abort(1);
	}
	
	if( QOP_F3_MDWF_import_fermion(&sprec_x0, 
				       state,
				       fermionReader,
				       (void *)&psi) != 0 ) {
	  QDPIO::cerr << "MDWF Error:  "<< QOP_MDWF_error(state) << endl;
	  QDP_abort(1);
	}
	
	if( QOP_F3_MDWF_allocate_fermion(&sprec_soln, 
					 state) != 0 ) {
	  QDPIO::cerr << "MDWF Error:  "<< QOP_MDWF_error(state) << endl;
	  QDP_abort(1);
	}
	
	/* Do the solve */
	
	double target_epsilon = toDouble(invParam.RsdCG*invParam.RsdCG);
	int max_iteration = invParam.MaxCG;
	int status; 
	
	QDPIO::cout << "MDWFQpropT: Beginning Single Precision Solve" << endl;
	if( ( status=QOP_F3_MDWF_DDW_CG(sprec_soln,
					&out_iters_single,
					&out_eps_single,
					params, 
					sprec_x0,
					sprec_gauge,
					sprec_rhs,
					max_iteration,
					target_epsilon, 
					QOP_MDWF_FINAL_DIRAC_RESIDUAL)) != 0 ) {
	  QDPIO::cerr << "MDWF Error:  "<< QOP_MDWF_error(state) << endl;
	  QDP_abort(1);
	}
	
	QDPIO::cout << "MDWF Single Prec : status =" << status 
		    << " iterations=" << out_iters_single
		    << " resulting epsilon=" << sqrt(out_eps_single) << endl;
	
	/* Export the solution */
	if( QOP_F3_MDWF_export_fermion(fermionWriter, 
				       &psi,
				       sprec_soln) != 0 ) { 
	  QDPIO::cerr << "MDWF Error:  "<< QOP_MDWF_error(state) << endl;
	  QDP_abort(1);
	}
	
	/* Now I can free the fermions */
	QOP_F3_MDWF_free_fermion(&sprec_soln);
	QOP_F3_MDWF_free_fermion(&sprec_x0);
	QOP_F3_MDWF_free_fermion(&sprec_rhs);
	QOP_F3_MDWF_free_gauge(&sprec_gauge);
	
	res.n_count = out_iters_single;

      }
#endif 
#if 1
      /* DoublePrecision Branch. Test this first
	 then duplicate for double prec */
      {
	int out_iters_double;
	double out_eps_double;

	// OK Try to solve for the single precision 
	// Cast to single precision gauge field
	QOP_D3_MDWF_Gauge *dprec_gauge = NULL;
	if ( QOP_D3_MDWF_import_gauge(&dprec_gauge,
				      state,
				      gaugeReader,
				      (void *)&u) != 0 ) { 
	  QDPIO::cerr << "MDWF Error:  "<< QOP_MDWF_error(state) << endl;
	  QDP_abort(1);
	}
	
	/* Fermion fields */
	QOP_D3_MDWF_Fermion *dprec_rhs;
	QOP_D3_MDWF_Fermion *dprec_x0;
	QOP_D3_MDWF_Fermion *dprec_soln;
	
	if( QOP_D3_MDWF_import_fermion(&dprec_rhs, 
				       state,
				       fermionReader,
				       (void *)&chi) != 0 ) {
	  QDPIO::cerr << "MDWF Error:  "<< QOP_MDWF_error(state) << endl;
	  QDP_abort(1);
	}
	
	if( QOP_D3_MDWF_import_fermion(&dprec_x0, 
				       state,
				       fermionReader,
				       (void *)&psi) != 0 ) {
	  QDPIO::cerr << "MDWF Error:  "<< QOP_MDWF_error(state) << endl;
	  QDP_abort(1);
	}
	
	if( QOP_D3_MDWF_allocate_fermion(&dprec_soln, 
					 state) != 0 ) {
	  QDPIO::cerr << "MDWF Error:  "<< QOP_MDWF_error(state) << endl;
	  QDP_abort(1);
	}
	
	/* Do the solve */
	
	double target_epsilon = toDouble(invParam.RsdCGRestart*invParam.RsdCGRestart);
	int max_iteration = invParam.MaxCGRestart;
	int status; 
	
	QDPIO::cout << "MDWFQpropT: Beginning Double Precision Solve" << endl;
	if( ( status=QOP_D3_MDWF_DDW_CG(dprec_soln,
					&out_iters_double,
					&out_eps_double,
					params, 
					dprec_x0,
					dprec_gauge,
					dprec_rhs,
					max_iteration,
					target_epsilon, 
					QOP_MDWF_LOG_EVERYTHING)) != 0 ) {
	  QDPIO::cerr << "MDWF Error:  "<< QOP_MDWF_error(state) << endl;
	  QDP_abort(1);
	}
	
	QDPIO::cout << "MDWF Double Prec : status =" << status 
		    << " iterations=" << out_iters_double
		    << " resulting epsilon=" << sqrt(out_eps_double) << endl;
	
	/* Export the solution */
	if( QOP_D3_MDWF_export_fermion(fermionWriter, 
				       &psi,
				       dprec_soln) != 0 ) { 
	  QDPIO::cerr << "MDWF Error:  "<< QOP_MDWF_error(state) << endl;
	  QDP_abort(1);
	}
	
	/* Now I can free the fermions */
	QOP_D3_MDWF_free_fermion(&dprec_soln);
	QOP_D3_MDWF_free_fermion(&dprec_x0);
	QOP_D3_MDWF_free_fermion(&dprec_rhs);
	QOP_D3_MDWF_free_gauge(&dprec_gauge);

	res.n_count += out_iters_double;
      }


#endif
      
      // Compute actual residual
      {
	multi1d<LatticeFermion>  r(N5);
	A->unprecLinOp(r, psi, PLUS);
	r -= chi;
	res.resid = sqrt(norm2(r));
      }
      QDPIO::cout << "MDWF Final: single_iters=" << out_iters_single << " double_iters=" << out_iters_double << " total_iters=" << res.n_count << endl;
      QDPIO::cout << "MDWF Final: final absolute unprec residuum="<<res.resid<<endl;

      return res;
    }


    /* INIT Function */
    void init(Handle< FermState<T,P,Q> > fermstate, const GroupXML_t& inv) 
    {
      if( Nd != 4 ) {
	QDPIO::cout << "This will only work for Nd=4" << endl;
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
      
      
      // I need to call Andrews init function.
      // For this I need a function that can tell me the 
      multi1d<int> lattice(5);
      multi1d<int> network(4);
      multi1d<int> node_coords(4);
      int master_p;
      
      // Lattice is the 5D lattice size
      for(int mu=0; mu < Nd; mu++) { 
	lattice[mu] = Layout::lattSize()[mu];
      }
      lattice[Nd]=N5;
      
      for(int mu=0; mu < Nd; mu++) { 
	network[mu] = Layout::logicalSize()[mu];
	node_coords[mu] = Layout::nodeCoord()[mu];
      }
      
      // Master_p has to be zero on master node and nonzero
      // elsewhere. This is odd.
      if( Layout::primaryNode() ) { 
	master_p = 0;
      }
      else { 
	master_p = 1;
      }
      
      // OK. Let's call Andrew's routine
      if(  QOP_MDWF_init(&state, lattice.slice(), network.slice(),
			 node_coords.slice(), master_p, sublattice_func, 
			 NULL) != 0 ) { 
	// Nonzero return value => error
	QDPIO::cerr << "MDWF Error: " << QOP_MDWF_error(state) << endl;
	QDP_abort(1);
      }
      
      // Set up the masses etc...
      u.resize(Nd);
      u = fermstate->getLinks();
      Real ff = Real(1);
      
      if (anisoParam.anisoP) {
	ff = where(anisoParam.anisoP, anisoParam.nu / anisoParam.xi_0, Real(1));
	for(int mu=0; mu < u.size(); ++mu) {
	  if (mu != anisoParam.t_dir)
	    u[mu] *= ff;
	}
      }
      
      // Set the Shamir parameters for now
      {
	double a5 = (double)1;

	// Convention change. Now in the internal code it is 5 - OverMass
	// To allow using anisotropy, I subtract off the 5 and 'add it back on'
	// suitable for the anisotropic case
	double M5  = -5 + toDouble(1 + a5*(1 + (Nd-1)*ff - OverMass));
	
	double m_f = toDouble(Mass);
	
	if(  QOP_MDWF_set_Shamir(&params, state, a5, M5 ,m_f) != 0){
	  QDPIO::cerr << "MDWF Error: " << QOP_MDWF_error(state)<< endl;
	  QDP_abort(1);
	}
	
      }

      return;
    }
      
    // Finalize - destructor call
    void fini(void)  {
      if (params != NULL) { 
	QOP_MDWF_free_parameters(&params);
      }
      
      

      if (state != NULL) { 
	QOP_MDWF_fini(&state);
      }

      
      return;
    }
    
    int size() const { return N5; }
    const Subset& subset() const { return all; }

    
  private:
    Handle< EvenOddPrecConstDetLinearOperatorArray<T,P,Q> > A;
    Real OverMass;              // M5
    Real Mass;                  // m_f
    int N5;                     // The 5th Dimension
    AnisoParam_t anisoParam;    // Anisotropy
    SysSolverCGParams invParam; // Inverter Parameters
    multi1d<LatticeColorMatrix> u; // The gauge field suitably prepared
    
    // Internal Pointers
    QOP_MDWF_State *state;
    QOP_MDWF_Parameters *params;
  };
  

  /*! @} */   // end of group qprop
}


#endif
