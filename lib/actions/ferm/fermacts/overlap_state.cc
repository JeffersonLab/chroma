#include "chromabase.h"
#include "actions/ferm/fermacts/overlap_state.h"

using namespace std;
using namespace QDP;


namespace Chroma {

  using namespace Chroma;
  using namespace OverlapConnectStateEnv;

  namespace OverlapConnectStateEnv {
 
    //! Create a ConnectState with just the gauge fields
    /*! Assumes that approximation bounds are those of the free field
        WilsonOp - ie. ApproxMin=0, ApproxMax=2Nd
    */
    template<typename T>
    static
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_, 
		       const FermBC<T>& fbc)
    {
      multi1d<LatticeColorMatrix> u_tmp = u_;
      fbc.modifyU(u_tmp);

      Real approxMin = 0.0;
      Real approxMax = Real(2)*Real(Nd);
      return new OverlapConnectState(u_tmp, approxMin, approxMax);
    }

    //! Create a ConnectState with just the gauge fields, and a lower
    //  approximation bound
    /*! Assumes that maximum approximation bounds are those of the free field
      WilsonOp - ie. ApproxMax=2Nd
    */
    template<typename T>
    static
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC<T>& fbc,
		       const Real& approxMin_)
    {
      if ( toBool( approxMin_ < Real(0) )) { 
	ostringstream error_str;
	error_str << "createOverlapState: approxMin_ has to be positive" 
		  << endl;
	throw error_str.str();
      }
 
      // First put in the BC
      multi1d<LatticeColorMatrix> u_tmp = u_;
      fbc.modifyU(u_tmp);
      
      // Set Approx Max
      Real approxMax = Real(2)*Real(Nd);
      return new OverlapConnectState(u_tmp, approxMin_, approxMax);
    }


    //! Create a connect State with just approximation range bounds
    template<typename T>
    static
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC<T>& fbc, 
		       const Real& approxMin_,
		       const Real& approxMax_)
    {
      ostringstream error_str;
      
      
      if ( toBool(approxMin_ < 0 )) { 
	error_str << "OverlapConnectState::createState: approxMin_ has to be positive" << endl;
	throw error_str.str();
      }
      
      if ( toBool(approxMax_ < approxMin_) ) { 
	error_str << "OverlapConnectState::createState: approxMax_ has to be larger than approxMin_" << endl;
	throw error_str.str();
      }
 
      
      // First put in the BC
      multi1d<LatticeColorMatrix> u_tmp = u_;
      fbc.modifyU(u_tmp);

      return new OverlapConnectState(u_tmp, approxMin_, approxMax_);
    }

    //! Create OverlapConnectState with eigenvalues/vectors
    template<typename T>
    static
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		const FermBC<T>& fbc, 
		const multi1d<Real>& lambda_lo_, 
		const multi1d<LatticeFermion>& evecs_lo_,
		const Real& lambda_hi_)
    {
      ostringstream error_str;
      
      if ( lambda_lo_.size() == 0 ) {
	error_str << "Attempt to createState with 0 e-values and no approxMin" << endl;
	throw error_str.str();
      }
      
      if ( lambda_lo_.size() != evecs_lo_.size() ) {
	error_str << "Attempt to createState with no of low eigenvalues != no of low eigenvectors" << endl;
	throw error_str.str();
      }
      
      Real approxMax = lambda_hi_;
      Real approxMin = fabs(lambda_lo_[ lambda_lo_.size() - 1 ]);
      
      // First put in the BC
      multi1d<LatticeColorMatrix> u_tmp = u_;
      fbc.modifyU(u_tmp);
      
      return new OverlapConnectState(u_tmp, 
				     lambda_lo_, 
				     evecs_lo_, 
				     lambda_hi_, 
				     approxMin, 
				     approxMax);
    }    


  
    //! Create a ConnectState out of XML
    /*! Needs the auxiliary fermion action linear operator 
        passed just in case it needs to check eigenvalues/vectors.
	The caller should be a fermact so it should be able to 
	generate the linop easy enough */
    template<typename T>
    static 
    const OverlapConnectState*  
    createOverlapState(const multi1d<LatticeColorMatrix> u_,
		       const FermBC<T>& fbc,
		       XMLReader& state_info_xml, 
		       const string& state_info_path,
		       const LinearOperator<LatticeFermion>& H)
    {
      OverlapStateInfo tmp_info;
      read(state_info_xml, state_info_path, tmp_info);
      return OverlapConnectStateEnv::createOverlapState<T>(u_, fbc, tmp_info, H);
    }

    //! Create from OverlapStateInfo Structure
    /*! Needs the auxiliary fermion action linear operator 
        passed just in case it needs to check eigenvalues/vectors.
	The caller should be a fermact so it should be able to 
	generate the linop easy enough */
    template<typename T>
    static 
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		const FermBC<T>& fbc,
		const OverlapStateInfo& state_info,
		const LinearOperator<LatticeFermion>& H)
    {
      
      // If No eigen values specified use min and max
      if ( state_info.getNWilsVec() == 0 ) { 
	return OverlapConnectStateEnv::createOverlapState<T>(u_,
						      fbc,
						      state_info.getApproxMin(),
						      state_info.getApproxMax());
      }
      else {
	
	// If there are eigen values, either load them, 
	if( state_info.loadEigVec() ) {
	  ChromaWilsonRitz_t ritz_header;
	  multi1d<Real> lambda_lo;
	  multi1d<LatticeFermion> eigv_lo;
	  Real lambda_hi;
	  const EigenIO_t& eigen_io = state_info.getEigenIO();
	  
	  if( eigen_io.eigen_filefmt == EVEC_TYPE_SCIDAC ) { 
	    readEigen(ritz_header, lambda_lo, eigv_lo, lambda_hi, 
		      eigen_io.eigen_file,
		      state_info.getNWilsVec(),
		      QDPIO_SERIAL);
	  }
	  else {
	    QDPIO::cerr << "Unsupported Eigenvector format for reading " << endl;
	    QDP_abort(1);
	  }
	  
	  QDPIO::cout << "createOverlapState: " << endl;
	  for(int i=0; i < lambda_lo.size(); i++) { 
	    QDPIO::cout << "lambda_lo["<<i<<"]= " << lambda_lo[i] << endl;
	  }
	  
	  QDPIO::cout << "|lambda_high|= " << lambda_hi << endl;
	  
	  // Test the e-values
	  // BEASTLY HACKERY!!!!
	  //  In order to test the evecs I need to create a ConnectState
	  //  for the fermions. I am assuming here, that the AuxiliaryFermAct
	  //  needs only a SimpleConnectState and I manufacture it by 
	  //  hand after applying the BC's of the calling Operator.
	  //  This goes hand in hand with the problem of turning 
	  //  a potentially 5D FermBC into a 4D one. If I could do that
	  //  then I wouldn't need to hack trivial 

	  multi1d<LatticeColorMatrix> u_test = u_;
	  fbc.modifyU(u_test);
	  Handle< const ConnectState > wils_connect_state(new SimpleConnectState(u_test));
	  
	  multi1d<Double> check_norm(state_info.getNWilsVec());
	  multi1d<Double> check_norm_rel(state_info.getNWilsVec());
	  for(int i=0; i < state_info.getNWilsVec() ; i++) { 
	    LatticeFermion Me;
	    H(Me, eigv_lo[i], PLUS);
	    
	    LatticeFermion lambda_e;
	    
	    lambda_e = lambda_lo[i]*eigv_lo[i];
	    LatticeFermion r_norm = Me - lambda_e;
	    check_norm[i] = sqrt(norm2(r_norm));
	    check_norm_rel[i] = check_norm[i]/fabs(Double(lambda_lo[i]));
	    
	    QDPIO::cout << "Eigenpair " << i << " Resid Norm = " 
			<< check_norm[i] << " Resid Rel Norm = " << check_norm_rel[i] << endl;
	  }

	  return OverlapConnectStateEnv::createOverlapState<T>(u_, fbc, lambda_lo, eigv_lo, lambda_hi);
	}
	else if( state_info.computeEigVec() ) {
	  QDPIO::cerr << "Recomputation of eigensystem not yet implemented" << endl;
	  QDP_abort(1);
	}
	else {
	  QDPIO::cerr << "I have to create a state without min/max, loading/computing eigenvectors/values. How do I do that? "<< endl;
	  QDP_abort(1);
	}
      }
      
      return 0;
    }


    // Non templated wrappers
    //! Create a ConnectState with just the gauge fields
    /*! Assumes that approximation bounds are those of the free field
        WilsonOp - ie. ApproxMin=0, ApproxMax=2Nd
    */
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_, 
		       const FermBC<LatticeFermion>& fbc) 
    {
      return OverlapConnectStateEnv::createOverlapState<LatticeFermion>(u_, fbc);
    }

    //! Create a ConnectState with just the gauge fields
    /*! Assumes that approximation bounds are those of the free field
        WilsonOp - ie. ApproxMin=0, ApproxMax=2Nd
    */
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_, 
		       const FermBC< multi1d<LatticeFermion> >& fbc) 
    {
      return OverlapConnectStateEnv::createOverlapState< multi1d<LatticeFermion> >(u_, fbc);
    }

    //! Create a ConnectState with just the gauge fields, and a lower
    //  approximation bound
    /*! Assumes that maximum approximation bounds are those of the free field
      WilsonOp - ie. ApproxMax=2Nd
    */
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC<LatticeFermion>& fbc,
		       const Real& approxMin_) 
    {
      return OverlapConnectStateEnv::createOverlapState<LatticeFermion>(u_, fbc, approxMin_ );
    }

    //! Create a ConnectState with just the gauge fields, and a lower
    //  approximation bound
    /*! Assumes that maximum approximation bounds are those of the free field
      WilsonOp - ie. ApproxMax=2Nd
    */
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC<multi1d<LatticeFermion> >& fbc,
		       const Real& approxMin_) 
    {
      return OverlapConnectStateEnv::createOverlapState< multi1d<LatticeFermion>  >(u_,fbc, approxMin_);
    }

    //! Create a connect State with just approximation range bounds
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC<LatticeFermion>& fbc, 
		       const Real& approxMin_,
		       const Real& approxMax_)
    {
      return OverlapConnectStateEnv::createOverlapState< LatticeFermion >(u_,fbc,approxMin_, approxMax_);
    }

    //! Create a connect State with just approximation range bounds
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC< multi1d<LatticeFermion> >& fbc, 
		       const Real& approxMin_,
		       const Real& approxMax_) 
    {
      return OverlapConnectStateEnv::createOverlapState<  multi1d<LatticeFermion> >(u_,fbc,approxMin_, approxMax_);
    }

    //! Create OverlapConnectState with eigenvalues/vectors
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC<LatticeFermion>& fbc, 
		       const multi1d<Real>& lambda_lo_, 
		       const multi1d<LatticeFermion>& evecs_lo_,
		       const Real& lambda_hi_) 
    {
      return OverlapConnectStateEnv::createOverlapState< LatticeFermion >(u_, fbc, lambda_lo_, evecs_lo_, lambda_hi_);
    }

    //! Create OverlapConnectState with eigenvalues/vectors
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC<multi1d<LatticeFermion> >& fbc, 
		       const multi1d<Real>& lambda_lo_, 
		       const multi1d<LatticeFermion>& evecs_lo_,
		       const Real& lambda_hi_) 
    {
      return OverlapConnectStateEnv::createOverlapState< multi1d<LatticeFermion>  >(u_, fbc, lambda_lo_, evecs_lo_, lambda_hi_);
    }


    //! Create a ConnectState out of XML
    /*! Needs the auxiliary fermion action linear operator 
        passed just in case it needs to check eigenvalues/vectors.
	The caller should be a fermact so it should be able to 
	generate the linop easy enough */
    const OverlapConnectState*  
    createOverlapState(const multi1d<LatticeColorMatrix> u_,
		       const FermBC< LatticeFermion >& fbc,
		       XMLReader& state_info_xml, 
		       const string& state_info_path,
		       const LinearOperator<LatticeFermion>& H)
    {
      return OverlapConnectStateEnv::createOverlapState< LatticeFermion >(u_, fbc, state_info_xml, 
						  state_info_path, H);
    }


    //! Create a ConnectState out of XML
    /*! Needs the auxiliary fermion action linear operator 
        passed just in case it needs to check eigenvalues/vectors.
	The caller should be a fermact so it should be able to 
	generate the linop easy enough */
    const OverlapConnectState*  
    createOverlapState(const multi1d<LatticeColorMatrix> u_,
		       const FermBC< multi1d< LatticeFermion > >& fbc,
		       XMLReader& state_info_xml, 
		       const string& state_info_path,
		       const LinearOperator<LatticeFermion>& H)
    {
      return OverlapConnectStateEnv::createOverlapState< multi1d< LatticeFermion > >(u_, fbc, state_info_xml, state_info_path, H);
    }


    //! Create from OverlapStateInfo Structure
    /*! Needs the auxiliary fermion action linear operator 
        passed just in case it needs to check eigenvalues/vectors.
	The caller should be a fermact so it should be able to 
	generate the linop easy enough */
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC< LatticeFermion >& fbc,
		       const OverlapStateInfo& state_info,
		       const LinearOperator<LatticeFermion>& H) 
    {
      return OverlapConnectStateEnv::createOverlapState< LatticeFermion >(u_, fbc, state_info, H);
    }

    //! Create from OverlapStateInfo Structure
    /*! Needs the auxiliary fermion action linear operator 
        passed just in case it needs to check eigenvalues/vectors.
	The caller should be a fermact so it should be able to 
	generate the linop easy enough */
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC< multi1d<LatticeFermion> >& fbc,
		       const OverlapStateInfo& state_info,
		       const LinearOperator<LatticeFermion>& H) 
    {
      return OverlapConnectStateEnv::createOverlapState< multi1d<LatticeFermion> >(u_, fbc, state_info, H);
    }

  }; // namespace OverlapConnectStateEnv
}; // namespace Chroma
