// -*- C++ -*-
// $Id: overlap_state.h,v 1.4 2004-09-27 14:58:43 bjoo Exp $
/*! @file
 * @brief Connection state holding eigenvectors
 *
 * Holds gauge fields and eigenvectors for overlap-ish thingies
 */

#ifndef __overlap_state_h__
#define __overlap_state_h__

#include "state.h"
#include "chromabase.h"
#include "io/overlap_state_info.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"

using namespace QDP;

namespace Chroma
{
  template<typename T>
  class OverlapConnectState
  {
  public:

    typedef Real WordBase_t;
  
    //! Constructor with no eigenvalues
    OverlapConnectState(const multi1d<LatticeColorMatrix>& u_,  // gauge field
			const WordBase_t& approxMin_,          // epsilon
			const WordBase_t& approxMax_           // approx max
      )  {
      u = u_;
      eigVal.resize(0);
      eigVec.resize(0);
      eigValMax = 0;
      approxMin = approxMin_ ;
      approxMax = approxMax_ ;
    }

    //! Constructor with e-values and e-vectors
    OverlapConnectState(const multi1d<LatticeColorMatrix>& u_,
			const multi1d<WordBase_t>& val_, 
			const multi1d<T>& vec_,
			const WordBase_t& val_max_,
			const WordBase_t& approxMin_,
			const WordBase_t& approxMax_)  {

      // hhmm, for now make a copy
      u = u_;
      eigVal = val_;
      eigVec = vec_;
      eigValMax = val_max_;
      approxMin = approxMin_;
      approxMax = approxMax_;
    }

    OverlapConnectState(const OverlapConnectState& a) : u(a.u), eigVal(a.eigVal), eigVec(a.eigVec), eigValMax(a.eigValMax), approxMin(a.approxMin), approxMax(a.approxMax) {}

    ~OverlapConnectState() {};

    //! Return the link fields needed in constructing linear operators
    const multi1d<LatticeColorMatrix>& getLinks() const {return u;}
  

    //! Return the eigenvalues
    const multi1d<WordBase_t>& getEigVal() const {return eigVal;}

    //! Return the eigenvectors
    const multi1d<T>& getEigVec() const {return eigVec;}
  
    //! Return the max eigenvalues
    const WordBase_t& getEigValMax() const {return eigValMax;}

    const WordBase_t& getApproxMin() const { return approxMin; }
    const WordBase_t& getApproxMax() const { return approxMax; }

  private:
    OverlapConnectState() {}  // hide default constructur
    void operator=(const OverlapConnectState&) {} // hide =

  private:
    multi1d<LatticeColorMatrix> u;
    multi1d<WordBase_t>  eigVal;
    multi1d<T>           eigVec;
    WordBase_t           eigValMax;
    WordBase_t           approxMin;
    WordBase_t           approxMax;

  };

  using namespace Chroma;
  namespace OverlapConnectStateEnv { 

    // Shared manipulation functions -- should only be called by fermacts

    //! Create a connect State with just approximation range bounds
    template <typename T> 
    static
    const OverlapConnectState<T>*
    createState(const multi1d<LatticeColorMatrix>& u_,
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

      return new OverlapConnectState<T>(u_tmp, approxMin_, approxMax_);
    }
    

    //! Create OverlapConnectState with eigenvalues/vectors
    template<typename T>
    static
    const OverlapConnectState<T>*
    createState(const multi1d<LatticeColorMatrix>& u_,
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
      
      return new OverlapConnectState<LatticeFermion>(u_tmp, lambda_lo_, evecs_lo_, lambda_hi_, approxMin, approxMax);
    }    


    //! Create from 
    template<typename T>
    static const OverlapConnectState<T>*
    createState(const multi1d<LatticeColorMatrix>& u_,
		const FermBC<T>& fbc,
		const OverlapStateInfo& state_info,
		const Real AuxMass) const
    {
      
      // If No eigen values specified use min and max
      if ( state_info.getNWilsVec() == 0 ) { 
	return createState(u_,
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
	  else if ( eigen_io.eigen_filefmt == EVEC_TYPE_SZIN ) { 
	    if( toBool( fabs(params.AuxMass) > 8 ) ){
	      QDPIO::cerr << "OverMass unspecified, or | OverMass | > 8" << endl;
	      QDPIO::cerr << "The wilson mass is needed to rescale the eigenvalues" << endl;
	      QDP_abort(1);
	    }
	    
	    readEigenSzin(lambda_lo, eigv_lo, lambda_hi, state_info.getNWilsVec(), eigen_io.eigen_file);
	    
	    // Now I need to scale things by the wilson mass (Nd + m)
	    for(int i=0; i < lambda_lo.size(); i++) { 
	      lambda_lo[i] *= (Real(Nd) + AuxMass);
	    }
	    
	    lambda_hi *= (Real(Nd) + AuxMass);
	    
	  }
	  else {
	    QDPIO::cerr << "Unsupported Eigenvector format for reading " << endl;
	    QDP_abort(1);
	  }
	  
	  QDPIO::cout << "createOverlapState: " << endl;
	  QDPIO::cout << " |lambda_lo|" << lambda_lo << endl;
	  QDPIO::cout << " |lambda_high|" << lambda_hi;
	  
	  // Test the e-values
	  multi1d<LatticeColorMatrix> u_test = u_;
	  fbc.modifyU(u_test);
	  Handle< const ConnectState > wils_connect_state(new SimpleConnectState(u_test));
	  Handle< const LinearOperator<LatticeFermion> > H = Mact->gamma5HermLinOp(wils_connect_state); 
	  
	  
	  multi1d<Double> check_norm(state_info.getNWilsVec());
	  multi1d<Double> check_norm_rel(state_info.getNWilsVec());
	  for(int i=0; i < state_info.getNWilsVec() ; i++) { 
	    LatticeFermion Me;
	    (*H)(Me, eigv_lo[i], PLUS);
	    
	    LatticeFermion lambda_e;
	    
	    lambda_e = lambda_lo[i]*eigv_lo[i];
	    LatticeFermion r_norm = Me - lambda_e;
	    check_norm[i] = sqrt(norm2(r_norm));
	    check_norm_rel[i] = check_norm[i]/fabs(Double(lambda_lo[i]));
	    
	    QDPIO::cout << "Eigenpair " << i << " Resid Norm = " 
			<< check_norm[i] << " Resid Rel Norm = " << check_norm_rel[i] << endl;
	  }
	  /*
	    write(xml_out, "eigen_norm", check_norm);
	    write(xml_out, "eigen_rel_norm", check_norm_rel);
	  */
	  
	  return createState(u_, fbc, lambda_lo, eigv_lo, lambda_hi);
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
  
    
    template<typename T>
    static 
    OverlapConnectState<T>*  
    createOverlapState(const multi1d<LatticeColorMatrix> u_,
		       const FermBC<T>& fbc,
		       XMLReader& state_info_xml, 
		       const string& state_info_path,
		       const Real AuxMass)
    {
      ZolotarevStateInfo tmp_info;
      read(state_info_xml, state_info_path, tmp_info);
      return OverlapConnectStateEnv::createState(u_, fbc, tmp_info, AuxMass);
    }
    
  };
  
};

#endif
