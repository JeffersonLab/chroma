// $Id: syssolver_linop_qdp_mg.cc, v1.0 2013-06-20 22:12 sdcohen $
/*! \file
 *  \brief Make contact with the QDP clover multigrid solver, transfer
 *         the gauge field, generate the coarse grids, solve systems
 */
#include "state.h"
#include "meas/inline/io/named_objmap.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include <cstdio>
#include <ostream>

using namespace std;

#if BASE_PRECISION == 32
#define QDP_Precision 'F'
#define QLA_Precision 'F'
#define toReal toFloat
#elif BASE_PRECISION == 64
#define QDP_Precision 'D'
#define QLA_Precision 'D'
#define toReal toDouble
#endif

#include "actions/ferm/invert/qop_mg/syssolver_linop_qop_mg_w.h"

extern "C" {
  // This should be placed on the include path.
#include "wilsonmg-interface.h"
extern struct MGP(Clover_Params) PC(g_param);
}

#include "meas/glue/mesplq.h"

namespace Chroma
{



  static multi1d<LatticeColorMatrix> u;
// These functions will allow QDP to look into the Chroma gauge field and set
// the QDP gauge field at each site equal to the one in Chroma. There doesn't
// seem to be a good way to treat the extra std::vector index of the gauge field,
// so we make a std::vector of functions to carry that index.
#define index(c) c[0]+nrow[0]/2*(c[1]+nrow[1]*(c[2]+nrow[2]*(c[3]+nrow[3]*((c[0]+c[1]+c[2]+c[3])%2)))) // Hrm, this didn't seem to work; using linearSiteIndex for now...
#define pepo(d) \
  void peekpoke##d(QLA(ColorMatrix) *dest, int coords[]) {\
    START_CODE();\
    multi1d<int> x(4); for (int i=0; i<4; i++) x[i] = coords[i];\
    ColorMatrix U; U.elem() = u[d].elem(Layout::linearSiteIndex(x));\
    END_CODE();\
    \
    QLA(Complex) z;\
    QLA(Real) real, imag;\
    for (int c1=0; c1<3; c1++)\
      for (int c2=0; c2<3; c2++) {\
        real = U.elem().elem().elem(c1,c2).real();\
        imag = U.elem().elem().elem(c1,c2).imag();\
        if (0&&c1==0&&c2==0) {std::fflush(stdout);printf("Chroma: gauge[%d] I am node %d, parsing %d %d %d %d; I see %g + I %g\n",d,Layout::nodeNumber(),x[0],x[1],x[2],x[3],real,imag); std::fflush(stdout);}\
        QLA(C_eq_R_plus_i_R)(&z, &real, &imag);\
        QLA_elem_M(*dest,c1,c2) = z;\
      }\
  }
  pepo(0)
  pepo(1)
  pepo(2)
  pepo(3)
#undef pepo
#undef index

  LinOpSysSolverQOPMG::
    LinOpSysSolverQOPMG(Handle< LinearOperator<T> > A_,
			Handle< FermState<T,Q,Q> > state_, 
                        const SysSolverQOPMGParams& invParam_) : 
  A(A_), state(state_), invParam(invParam_), subspace(0x0)
  {
  // Copy the parameters read from XML into the QDP global structure
    for (int d=0; d<4; d++) PC(g_param).bc[d]  = 1;
    PC(g_param).aniso_xi = toReal(invParam.AnisoXi);
    PC(g_param).aniso_nu = toReal(invParam.AnisoNu);
    PC(g_param).kappa    = toReal(invParam.Kappa);
    PC(g_param).kappac   = toReal(invParam.KappaCrit);
    PC(g_param).mass     = toReal(invParam.Mass);
    PC(g_param).massc    = toReal(invParam.MassCrit);
    PC(g_param).clov_s   = toReal(invParam.Clover);
    PC(g_param).clov_t   = toReal(invParam.CloverT);
    PC(g_param).res      = toReal(invParam.Residual);
    PC(g_param).ngcr     = invParam.NumGCRVecs;
    PC(g_param).maxiter  = invParam.MaxIter;
    PC(g_param).verb     = invParam.Verbose;
    PC(g_param).levels   = invParam.Levels;

    for (int n=0; n<invParam.Levels; n++) {
      for (int d=0; d<4; d++)
        PC(g_param).block[n][d] = invParam.Blocking[n][d];
      PC(g_param).nNullVecs[n]   = invParam.NumNullVecs[n];
      PC(g_param).nullMaxIter[n] = invParam.NullMaxIter[n];
      PC(g_param).nullRes[n]  = toReal(invParam.NullResidual[n]);
      PC(g_param).nullConv[n] = toReal(invParam.NullConvergence[n]);
      PC(g_param).nExtraVecs[n]  = invParam.NumExtraVecs[n];
      PC(g_param).urelax[n]   = toReal(invParam.Underrelax[n]);
      PC(g_param).npre[n]        = invParam.NumPreHits[n];
      PC(g_param).npost[n]       = invParam.NumPostHits[n];
      PC(g_param).cngcr[n]       = invParam.CoarseNumGCRVecs[n];
      PC(g_param).cmaxiter[n]    = invParam.CoarseMaxIter[n];
      PC(g_param).cres[n]     = toReal(invParam.CoarseResidual[n]);
    }

// The std::vector of functions will be used by QDP to assign the gauge links
// at each site of the QDP lattice
    if (invParam.Levels>0) {
      /* We're going to pull the gauge field out of Chroma's aether */
#if 0
      u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(invParam.GaugeID);
#else
      u = state_->getLinks();
#endif
      // Compute the plaquette for comparison with MG code
      {
      	Double w_plaq, s_plaq, t_plaq, link;
	
      	MesPlq(u, w_plaq, s_plaq, t_plaq, link);
      	QDPIO::cout << "Plaquette from State: " << std::endl;
      	QDPIO::cout << "  w_plaq = " << w_plaq << std::endl;
      	QDPIO::cout << "  s_plaq = " << s_plaq << std::endl;
      	QDPIO::cout << "  t_plaq = " << t_plaq << std::endl;
      	QDPIO::cout << "  link trace =  " << link << std::endl;
      }
      
      int machsize[4], latsize[4];
      for (int d=0;d<4;d++) machsize[d] = Layout::logicalSize()[d];
      for (int d=0;d<4;d++) latsize[d]  = Layout::lattSize()[d];
      void (*peekpoke[4])(QLA(ColorMatrix) *dest, int coords[]) =
        {peekpoke0,peekpoke1,peekpoke2,peekpoke3};
      MGP(initialize)(machsize, latsize, peekpoke);
	
	//MGP(teststuff)();
    }
  }      

  LinOpSysSolverQOPMG::
    ~LinOpSysSolverQOPMG()
  {
    // If we don't use an Externa subspace we should free it here.
    if ( !invParam.ExternalSubspace ) { 
      if ( subspace != 0x0 ) MGP(destroy_subspace)(subspace);
    }
    MGP(finalize)();
  }

  
  void *fermionsrc, *fermionsol;
  template<typename T>
  void peekpokesrc(QLA(DiracFermion) *dest, int coords[])
  {
    multi1d<int> x(4); for (int i=0; i<4; i++) x[i] = coords[i];
    Fermion src; src.elem() = ((T*)fermionsrc)->elem(Layout::linearSiteIndex(x));
    /*START_CODE();
    double bsq = norm2(*(T*)fermionsrc).elem().elem().elem().elem();
    printf("Chroma:   in norm2 = %g\n",bsq);
    END_CODE();*/
    //printf("Chroma: x = %i %i %i %i:\n",x[0],x[1],x[2],x[3]);
    QLA(Complex) z;
    QLA(Real) real, imag;
    for (int s=0; s<4; s++)
      for (int c=0; c<3; c++) {
        real = src.elem().elem(s).elem(c).real();
        imag = src.elem().elem(s).elem(c).imag();
        //printf("Chroma:   s=%i,c=%i == %g + I %g\n",s,c,real,imag);
        QLA(C_eq_R_plus_i_R)(&z, &real, &imag);
        QLA(elem_D)(*dest,c,s) = z;
      }
  }
  template<typename T>
  void peekpokesol(QLA(DiracFermion) *src, int coords[])
  {
    multi1d<int> x(4); for (int i=0; i<4; i++) x[i] = coords[i];
    ColorVector ctmp;
    Fermion ftmp;
    /*START_CODE();
    double bsq = norm2(*(T*)fermionsrc).elem().elem().elem().elem();
    printf("Chroma:   in norm2 = %g\n",bsq);
    END_CODE();*/
    //printf("Chroma: x = %i %i %i %i:\n",x[0],x[1],x[2],x[3]);
    QLA(Complex) z;
    QLA(Real) real, imag;
    for (int s=0; s<4; s++) {
      for (int c=0; c<3; c++) {
        z = QLA_elem_D(*src,c,s);
        QLA(R_eq_re_C)(&real, &z);
        QLA(R_eq_im_C)(&imag, &z);
        Complex ztmp = cmplx(Real(real), Real(imag));
        pokeColor(ctmp, ztmp, c);
      }
      pokeSpin(ftmp, ctmp, s);
    }
    pokeSite(*((T*)fermionsol), ftmp, x);
  }

  //! Solve the linear system
  /*!
   * \param psi      solution ( Modify )
   * \param chi      source ( Read )
   * \return syssolver results
   */
  SystemSolverResults_t
    LinOpSysSolverQOPMG::operator() (T& psi, const T& chi) const
  {
    START_CODE();
    
    SystemSolverResults_t res;
    
    StopWatch swatch;
    swatch.reset();
    swatch.start();

    int machsize[4], latsize[4];
    for (int d=0;d<4;d++) machsize[d] = Layout::logicalSize()[d];
    for (int d=0;d<4;d++) latsize[d]  = Layout::lattSize()[d];

    // This will get the subspace
    //   If the subspace is externally managed it will be loaded/created
    //   If the subspace is internally managed it will be used/created
    WilsonMGSubspace* solve_subspace = getSubspace();


    // Set global pointers to our source and solution fermion fields
    fermionsrc = (void*)&chi;
    fermionsol = (void*)&psi;
    Double bsq = norm2(chi,all);
    QDPIO::cout << "Chroma:   chi all norm2 = " << bsq << std::endl;

    // DOING SOLVE HERE
    res.n_count = MGP(solve)(peekpokesrc<T>, peekpokesol<T>,solve_subspace);

    bsq = norm2(psi,all);
    QDPIO::cout << "Chroma:   psi all norm2 = " << bsq << std::endl;
    swatch.stop();
    Double rel_resid;
    double time = swatch.getTimeInSeconds();
    { 
      T r;
      r[A->subset()] = chi;
      T tmp;
      (*A)(tmp, psi, PLUS);
      r[A->subset()] -= tmp;
      res.resid = sqrt(norm2(r, A->subset()));
      rel_resid = res.resid / sqrt(norm2(chi, A->subset()));
    }

    QDPIO::cout << "QOPMG_LINOP_SOLVER: " 
		<< " Mass: " << invParam.Mass 
		<< " iterations: " << res.n_count
                << " time: " << time << " secs"
                << " Rsd: " << res.resid
                << " Relative Rsd: " << rel_resid
                << std::endl;

    if ( invParam.TerminateOnFail ) { 
      if ( toBool( rel_resid >= invParam.RsdToleranceFactor*invParam.Residual ) ) {
	QDPIO::cout << "!!!!! QOPMG_LINOP_SOLVER: " 
		    << "Mass: "  << invParam.Mass
		    << " did NOT CONVERGE: Target RelRsd = " << invParam.Residual
		    << " Actual RelRsd = " << rel_resid << std::endl;
	QDPIO::cout << "Exiting !!!!! " << std::endl;
	QDP_abort(1);
      }
    }
	
     
    END_CODE();
    return res;
  }

  LinOpSysSolverQOPMG::WilsonMGSubspace* LinOpSysSolverQOPMG::getSubspace() const
  {
    WilsonMGSubspace *ret_val=0x0;
    int machsize[4], latsize[4];
    for (int d=0;d<4;d++) machsize[d] = Layout::logicalSize()[d];
    for (int d=0;d<4;d++) latsize[d]  = Layout::lattSize()[d];

    QDPIO::cout << "LinOpSysSolverQOPMG::getSubspace() " << std::endl;
    if( !invParam.ExternalSubspace ) {
      if( !subspace ) { 
	QDPIO::cout << "Internal doesn't yet exist... creating" << std::endl;
	subspace = MGP(create_subspace)(latsize);
      }
      ret_val=subspace;
    }
    else {
      // Check if ID exists in the map.
      if( TheNamedObjMap::Instance().check(invParam.SubspaceId)  ) {
	QDPIO::cout << "Retreiving Subspace Pointer from NamedObject Map" << endl;	
	ret_val = TheNamedObjMap::Instance().getData< WilsonMGSubspace * >(invParam.SubspaceId);
	MGP(reset_subspace)(latsize,ret_val);
      }
      else { 
	QDPIO::cout << "Using External Subspace but object ID: " << invParam.SubspaceId << " is not present in the NamedObject store." << std:: endl;
	QDPIO::cout << "Creating a new one" << std::endl;
	ret_val = MGP(create_subspace)(latsize);
	
	QDPIO::cout << "Saving in Map" << std::endl;
	QDPIO::cout << "Creating Named Object Map Entry for subspace" << std::endl;
	XMLBufferWriter file_xml;
	push(file_xml, "FileXML");
	pop(file_xml);

	XMLBufferWriter record_xml;
	push(record_xml, "RecordXML");
	write(record_xml, "InvertParam", invParam);
	pop(record_xml);

	TheNamedObjMap::Instance().create< WilsonMGSubspace *>(invParam.SubspaceId);
	TheNamedObjMap::Instance().get(invParam.SubspaceId).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(invParam.SubspaceId).setRecordXML(record_xml);
	TheNamedObjMap::Instance().getData<WilsonMGSubspace*>(invParam.SubspaceId)=ret_val;
      } // else of if Subspace ID is in the map
    } // else of if !External Subspace

    return ret_val; // Return whichever pointer was the good one.
  }

  void  LinOpSysSolverQOPMG::eraseSubspace()  
  {
    QDPIO::cout << "LinOpSysSolverQOPMG: Erasing Subspace" << endl;
    if( !invParam.ExternalSubspace ) { 

      QDPIO::cout << "My subspace does not come from the NamedObjMap" << std::endl;

      // Internal subspace case, if subspace is not null then free it.
      if ( subspace != 0x0 ) {
	QDPIO::cout << "My subspace pointer is not null... Destroying" 
		    << std::endl;

	MGP(destroy_subspace)(subspace);
	subspace = 0x0;
      }
    }
    else { 
      // Using external subspace. If it exists free it.
      if ( TheNamedObjMap::Instance().check(invParam.SubspaceId) ) {
	QDPIO::cout << "MG Subspace: " << invParam.SubspaceId << " found. " << std::endl;
	
	// The loaded subspace ID is in the Map. Get it.
	WilsonMGSubspace *named_obj_subspace = TheNamedObjMap::Instance().getData<WilsonMGSubspace*>(invParam.SubspaceId);
	MGP(destroy_subspace)(named_obj_subspace);
	try { 
	  QDPIO::cout << "Attempting to delete from named object store" << std:: endl;
	  // Now erase the object
	  TheNamedObjMap::Instance().erase(invParam.SubspaceId);
	  
	  QDPIO::cout << "Object erased" << std::endl;
	}
	catch( std::bad_cast ) {
	  QDPIO::cerr <<" MGMdagMInternal::eraseSubspace: cast error" 
		      << std::endl;
	  QDP_abort(1);
	}
	catch (const std::string& e) {
	  
	  QDPIO::cerr << " MGMdagMInternal::eraseSubspace: error message: " << e 
		      << std::endl;
	  QDP_abort(1);
	}
	
	QDPIO::cout << "Subspace: " << invParam.SubspaceId << " destroyed and removed from map" << std::endl;
      }
      else { 
	QDPIO::cout << "MG Subspace: " << invParam.SubspaceId << " is not in the map. Nothing to desroy" << std::endl;
      }
    }
  }

  //! QDP multigrid system solver namespace
  namespace LinOpSysSolverQOPMGEnv
  {
    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("QOP_CLOVER_MULTIGRID_INVERTER");

      //! Local registration flag
      bool registered = false;
    }


    //! Callback function for standard precision
    LinOpSystemSolver<LatticeFermion>*
      createFerm( XMLReader& xml_in,
                  const std::string& path,
                  Handle< FermState< LatticeFermion, 
                                     multi1d<LatticeColorMatrix>,
                                     multi1d<LatticeColorMatrix> > > state, 
                  Handle< LinearOperator<LatticeFermion> >           A)
    {
      return new LinOpSysSolverQOPMG(A, state, SysSolverQOPMGParams(xml_in, path));
    }

    /*//! Callback function for single precision
    LinOpSystemSolver<LatticeFermionF>*
      createFermF( XMLReader&                                          xml_in,
                   const std::string&                                  path,
                   Handle< FermState< LatticeFermionF, 
                                      multi1d<LatticeColorMatrixF>,
                                      multi1d<LatticeColorMatrixF> > > state,
                   Handle< LinearOperator<LatticeFermionF> >           A)
    {
      return new LinOpSysSolverQOPMG<LatticeFermionF>
                   (A, SysSolverQOPMGParams(xml_in, path));
    }*/

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
        success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
        //success &= Chroma::TheLinOpFFermSystemSolverFactory::Instance().registerObject(name, createFermF);
        registered = true;
      }
      return success;
    }
  }
}
