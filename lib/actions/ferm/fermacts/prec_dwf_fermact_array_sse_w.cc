// $Id: prec_dwf_fermact_array_sse_w.cc,v 1.5 2004-10-18 20:40:57 edwards Exp $
/*! \file
 *  \brief SSE 4D style even-odd preconditioned domain-wall fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_array_sse_w.h"
#include "actions/ferm/linop/unprec_dwf_linop_array_w.h"
#include "actions/ferm/linop/prec_dwf_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include <sse_dwf_cg.h>

#include "actions/ferm/fermacts/fermfactory_w.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace SSEEvenOddPrecDWFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct< multi1d<LatticeFermion> >* createFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
								XMLReader& xml_in,
								const std::string& path)
    {
      return new SSEEvenOddPrecDWFermActArray(fbc, SSEEvenOddPrecDWFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    EvenOddPrecDWFermActBaseArray<LatticeFermion>* createDWFermAct(Handle< FermBC< multi1d<LatticeFermion> > > fbc,
								   XMLReader& xml_in,
								   const std::string& path)
    {
      return new SSEEvenOddPrecDWFermActArray(fbc, SSEEvenOddPrecDWFermActArrayParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "SSE_DWF";    // TEMPORARY HACK

    //! Register the Wilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActArrayFactory::Instance().registerObject(name, createFermAct)
                          & Chroma::TheEvenOddPrecDWFermActBaseArrayFactory::Instance().registerObject(name, createDWFermAct); 
  }


  //! Read parameters
  SSEEvenOddPrecDWFermActArrayParams::SSEEvenOddPrecDWFermActArrayParams(XMLReader& xml, 
									 const std::string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    read(paramtop, "OverMass", OverMass);
    read(paramtop, "Mass", Mass);
    read(paramtop, "N5", N5);

    if (paramtop.count("a5") != 0) 
      read(paramtop, "a5", a5);
    else
      a5 = 1.0;
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, SSEEvenOddPrecDWFermActArrayParams& param)
  {
    SSEEvenOddPrecDWFermActArrayParams tmp(xml, path);
    param = tmp;
  }



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
  const EvenOddPrecDWLinOpBaseArray<LatticeFermion>*
  SSEEvenOddPrecDWFermActArray::linOp(Handle<const ConnectState> state) const
  {
    return new EvenOddPrecDWLinOpArray(state->getLinks(),OverMass,Mass,N5);
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
  const UnprecDWLinOpBaseArray<LatticeFermion>*
  SSEEvenOddPrecDWFermActArray::linOpPV(Handle<const ConnectState> state) const
  {
    // For the PV operator, use the **unpreconditioned** one
    // fixed to quark mass 1
    return new UnprecDWLinOpArray(state->getLinks(),OverMass,1.0,N5);
  }


  //! Gauge field reader
  static double
  gauge_reader(const void *ptr, void *env, 
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
  static double
  fermion_reader(const void *ptr, void *env, 
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
  static double
  fermion_reader_zero(const void *ptr, void *env, 
		      const int latt_coord[5], int color, int spin, int reim)
  {
    return 0.0;
  }


  //! Fermion field reader
  static void
  fermion_writer(void *ptr, void *env, 
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



  ///////////////////////////////////////////////////////////////////////////////
  static void
  solve_cg5(multi1d<LatticeFermion> &solution,    // output
	    const multi1d<LatticeColorMatrix> &U, // input
	    double M_0,                           // input
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
      V[i] = shift(U[i], -1, i); // as viewed from the destination (sic)

    SSE_DWF_Gauge *g = SSE_DWF_load_gauge(&U, &V, NULL, gauge_reader);
    SSE_DWF_Fermion *eta = SSE_DWF_load_fermion(&rhs, NULL, fermion_reader);
    SSE_DWF_Fermion *X0 = SSE_DWF_load_fermion(&x0, NULL, fermion_reader);
    SSE_DWF_Fermion *res = SSE_DWF_allocate_fermion();

    double out_eps;
    int status = SSE_DWF_cg_solver(res, &out_eps, &out_iter,
				   g, M_0, m_f, X0, eta, 
				   rsd, max_iter);

    QDPIO::cout << "SSE DWF solver: status = " << status
		<< ", iterations = " << out_iter
		<< ", resulting epsilon = " << out_eps
		<< endl;

    SSE_DWF_save_fermion(&solution, NULL, fermion_writer, res);

    SSE_DWF_delete_fermion(res);
    SSE_DWF_delete_fermion(X0);
    SSE_DWF_delete_fermion(eta);
    SSE_DWF_delete_gauge(g);

    // SSE_DWF_fini();
  }


  //! Optimized inverter - this is temporary
  void 
  SSEEvenOddPrecDWFermActArray::opt_qpropT(multi1d<LatticeFermion>& psi, 
					   Handle<const ConnectState> state, 
					   const multi1d<LatticeFermion>& chi, 
					   const InvertParam_t& invParam,
					   int& ncg_had) const
  {
    QDPIO::cout << "entering SSEEvenOddPrecDWFermActArray::opt_qpropT" << endl;

    START_CODE();

    const multi1d<LatticeColorMatrix>& u = state->getLinks();

    multi1d<LatticeFermion> avp_chi(N5),tmp5(N5),phi(N5);

    //apply chroma operator
    Handle<const LinearOperator< multi1d<LatticeFermion> > > dwf_chroma(new UnprecDWLinOpArray(u,OverMass,Mass,N5));
    (*dwf_chroma)(phi, chi, PLUS);
    //phi = dwf_chroma * chi 
    // avp_chi = R*gamma_5 * phi 
    for(int s(0);s<N5;s++)
      avp_chi[N5-1-s] = Gamma(15)*phi[s] ;

    // Apply SSE inverter
    double M0 = toDouble(-2*(5.0-OverMass));
    double m_f = toDouble(Mass);
    double rsd = toDouble(invParam.RsdCG);
    int    max_iter = invParam.MaxCG;
    solve_cg5(tmp5, u, M0, m_f, avp_chi, avp_chi,
	      rsd, max_iter, ncg_had);

    // Need remaping to interface with andrew R\gamma_5
    for(int s(0);s<N5;s++) {
      psi[s] = -2 * (Gamma(15)*tmp5[N5-1-s]);
    }

    END_CODE();

    QDPIO::cout << "exiting SSEEvenOddPrecDWFermActArray::opt_qpropT" << endl;
  }

}
