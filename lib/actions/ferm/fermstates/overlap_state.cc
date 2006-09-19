// $Id: overlap_state.cc,v 1.1 2006-09-19 17:53:37 edwards Exp $
/*! @file
 * @brief Connection state holding eigenvectors
 *
 * Holds gauge fields and eigenvectors for overlap-ish thingies
 */

#include "chromabase.h"
#include "actions/ferm/fermstates/overlap_state.h"
#include "actions/ferm/fermstates/simple_fermstate.h"


namespace Chroma 
{

  // Save some typing
  typedef FermBC<LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >  FBC;

  //! Constructor with no eigenvalues
  OverlapConnectState::OverlapConnectState(Handle<FBC> fbc_,
					   const multi1d<LatticeColorMatrix>& u_,  // gauge field
					   const WordBase_t& approxMin_,          // epsilon
					   const WordBase_t& approxMax_)         // approx max
  {
    init(fbc_, u_, approxMin_, approxMax_);
  }
  
  //! Constructor with e-values and e-vectors
  OverlapConnectState::OverlapConnectState(Handle<FBC> fbc_,
					   const multi1d<LatticeColorMatrix>& u_,
					   const multi1d<WordBase_t>& val_, 
					   const multi1d<LatticeFermion>& vec_,
					   const WordBase_t& val_max_,
					   const WordBase_t& approxMin_,
					   const WordBase_t& approxMax_)
  {
    init(fbc_, u_, val_, vec_, val_max_, approxMin_, approxMax_);
  }


  //! Constructor with no eigenvalues
  void OverlapConnectState::init(Handle<FBC> fbc_,
				 const multi1d<LatticeColorMatrix>& u_,  // gauge field
				 const WordBase_t& approxMin_,          // epsilon
				 const WordBase_t& approxMax_)          // approx max
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
 
    // hhmm, for now make a copy
    fbc = fbc_;
    u   = u_;
    fbc->modify(u);  // apply BC

    eigVal.resize(0);
    eigVec.resize(0);
    eigValMax = 0;
    approxMin = approxMin_;
    approxMax = approxMax_;
  }
  
  //! Constructor with e-values and e-vectors
  void OverlapConnectState::init(Handle<FBC> fbc_,
				 const multi1d<LatticeColorMatrix>& u_,
				 const multi1d<WordBase_t>& val_, 
				 const multi1d<LatticeFermion>& vec_,
				 const WordBase_t& val_max_,
				 const WordBase_t& approxMin_,
				 const WordBase_t& approxMax_)
  {
    ostringstream error_str;
      
    if ( val_.size() == 0 ) {
      error_str << "Attempt to createState with 0 e-values and no approxMin" << endl;
      throw error_str.str();
    }
      
    if ( val_.size() != vec_.size() ) {
      error_str << "Attempt to createState with no of low eigenvalues != no of low eigenvectors" << endl;
      throw error_str.str();
    }
      
    // hhmm, for now make a copy
    fbc = fbc_;
    u   = u_;
    fbc->modify(u);  // apply BC

    eigVal = val_;
    eigVec = vec_;
    eigValMax = val_max_;
    approxMin = approxMin_;
    approxMax = approxMax_;
  }


  //! Create a ConnectState with just the gauge fields
  /*! Assumes that approximation bounds are those of the free field
    WilsonOp - ie. ApproxMin=0, ApproxMax=2Nd
  */
  OverlapConnectState::OverlapConnectState(Handle<FBC> fbc_,
					   const multi1d<LatticeColorMatrix>& u_)
  {
    Real approxMin_ = 0.0;
    Real approxMax_ = Real(2)*Real(Nd);
    init(fbc_, u_, approxMin_, approxMax_);
  }

  //! Create a ConnectState with just the gauge fields, and a lower
  //  approximation bound
  /*! Assumes that maximum approximation bounds are those of the free field
    WilsonOp - ie. ApproxMax=2Nd
  */
  OverlapConnectState::OverlapConnectState(Handle<FBC> fbc_,
					   const multi1d<LatticeColorMatrix>& u_,
					   const Real& approxMin_)
  {
    // Set Approx Max
    Real approxMax_ = Real(2)*Real(Nd);
    init(fbc_, u_, approxMin_, approxMax_);
  }


  //! Create OverlapConnectState with eigenvalues/vectors
  OverlapConnectState::OverlapConnectState(Handle<FBC> fbc_,
					   const multi1d<LatticeColorMatrix>& u_,
					   const multi1d<Real>& lambda_lo_, 
					   const multi1d<LatticeFermion>& evecs_lo_,
					   const Real& lambda_hi_) 
  {
    Real approxMax_ = lambda_hi_;
    Real approxMin_ = fabs(lambda_lo_[ lambda_lo_.size() - 1 ]);

    init(fbc, u_, lambda_lo_, evecs_lo_, lambda_hi_, approxMin_, approxMax_);
  }



#if 1
  //! Create a ConnectState out of XML
  /*! Needs the auxiliary fermion action linear operator 
    passed just in case it needs to check eigenvalues/vectors.
    The caller should be a fermact so it should be able to 
    generate the linop easy enough */
  OverlapConnectState::OverlapConnectState(Handle<FBC> fbc_,
					   const multi1d<LatticeColorMatrix> u_,
					   XMLReader& state_info_xml, 
					   const string& state_info_path,
					   const LinearOperator<LatticeFermion>& H)
  {
    OverlapStateInfo tmp_info;
    read(state_info_xml, state_info_path, tmp_info);
    init(fbc_, u_, tmp_info, H);
  }
#endif


  //! Create from OverlapStateInfo Structure
  /*! Needs the auxiliary fermion action linear operator 
    passed just in case it needs to check eigenvalues/vectors.
    The caller should be a fermact so it should be able to 
    generate the linop easy enough */
  OverlapConnectState::OverlapConnectState(Handle<FBC> fbc_,
					   const multi1d<LatticeColorMatrix>& u_,
					   const OverlapStateInfo& state_info,
					   const LinearOperator<LatticeFermion>& H) 
  {
    init(fbc_, u_, state_info, H);
  }


  //! Create from OverlapStateInfo Structure
  /*! Needs the auxiliary fermion action linear operator 
    passed just in case it needs to check eigenvalues/vectors.
    The caller should be a fermact so it should be able to 
    generate the linop easy enough */
  void OverlapConnectState::init(Handle<FBC> fbc_,
				 const multi1d<LatticeColorMatrix>& u_,
				 const OverlapStateInfo& state_info,
				 const LinearOperator<LatticeFermion>& H) 
  {
    // If No eigen values specified use min and max
    if ( state_info.getNWilsVec() == 0 ) 
    { 
      init(fbc_, u_, state_info.getApproxMin(), state_info.getApproxMax());
    }
    else 
    {
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
	//  needs only a SimpleFermState and I manufacture it by 
	//  hand after applying the BC's of the calling Operator.
	//  This goes hand in hand with the problem of turning 
	//  a potentially 5D FermBC into a 4D one. If I could do that
	//  then I wouldn't need to hack trivial 
	//
	//  RGE: maybe this problem is solved. There is really no difference in 4D and 5D
	//  fermbcs anymore. Anyway, the fbc is passed into the SimpleFermState now.

	Handle< FermState<T,P,Q> > wils_connect_state(new SimpleFermState<T,P,Q>(fbc_, u_));
	  
	multi1d<Double> check_norm(state_info.getNWilsVec());
	multi1d<Double> check_norm_rel(state_info.getNWilsVec());
	for(int i=0; i < state_info.getNWilsVec() ; i++) 
	{ 
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

	Real approxMax_ = lambda_hi;
	Real approxMin_ = fabs(lambda_lo[ lambda_lo.size() - 1 ]);

	init(fbc_, u_, lambda_lo, eigv_lo, lambda_hi, approxMin_, approxMax_);
      }
      else if( state_info.computeEigVec() ) 
      {
	QDPIO::cerr << "Recomputation of eigensystem not yet implemented" << endl;
	QDP_abort(1);
      }
      else 
      {
	QDPIO::cerr << "I have to create a state without min/max, loading/computing eigenvectors/values. How do I do that? "<< endl;
	QDP_abort(1);
      }
    }
  }


} // namespace Chroma

