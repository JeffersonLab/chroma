// -*- C++ -*-
// $Id: overlap_state.h,v 1.7 2005-01-14 20:13:04 edwards Exp $
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
#include "fermact.h"



namespace Chroma
{
  class OverlapConnectState : public ConnectState
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
			const multi1d<LatticeFermion>& vec_,
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
    const multi1d<LatticeFermion>& getEigVec() const {return eigVec;}
  
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
    multi1d<LatticeFermion> eigVec;
    WordBase_t           eigValMax;
    WordBase_t           approxMin; 
    WordBase_t           approxMax;

  };

  using namespace Chroma;
  namespace OverlapConnectStateEnv { 
    // Shared manipulation functions -- should only be called by fermacts
    

    //! Create a ConnectState with just the gauge fields
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_, 
		       const FermBC<LatticeFermion>& fbc);

    //! Create a ConnectState with just the gauge fields, and a lower
    //!  approximation bound
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC<LatticeFermion>& fbc,
		       const Real& approxMin_);

    //! Create a connect State with just approximation range bounds
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC<LatticeFermion>& fbc, 
		       const Real& approxMin_,
		       const Real& approxMax_);

    //! Create OverlapConnectState with eigenvalues/vectors
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC<LatticeFermion>& fbc, 
		       const multi1d<Real>& lambda_lo_, 
		       const multi1d<LatticeFermion>& evecs_lo_,
		       const Real& lambda_hi_);

    //! Create a ConnectState out of XML
    const OverlapConnectState*  
    createOverlapState(const multi1d<LatticeColorMatrix> u_,
		       const FermBC< LatticeFermion >& fbc,
		       XMLReader& state_info_xml, 
		       const string& state_info_path,
		       const LinearOperator<LatticeFermion>& H);

    
    //! Create from OverlapStateInfo Structure
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC< LatticeFermion >& fbc,
		       const OverlapStateInfo& state_info,
		       const LinearOperator<LatticeFermion>& H);
    


    //! Create a ConnectState with just the gauge fields
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_, 
		       const FermBC< multi1d<LatticeFermion> >& fbc);

    //! Create a ConnectState with just the gauge fields, and a lower
    //!  approximation bound
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC< multi1d<LatticeFermion> >& fbc,
		       const Real& approxMin_);


    //! Create a connect State with just approximation range bounds
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC< multi1d<LatticeFermion> >& fbc, 
		       const Real& approxMin_,
		       const Real& approxMax_);


    //! Create OverlapConnectState with eigenvalues/vectors
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC< multi1d<LatticeFermion> >& fbc, 
		       const multi1d<Real>& lambda_lo_, 
		       const multi1d<LatticeFermion>& evecs_lo_,
		       const Real& lambda_hi_);

    //! Create a ConnectState out of XML
    const OverlapConnectState*  
    createOverlapState(const multi1d<LatticeColorMatrix> u_,
		       const FermBC< multi1d<LatticeFermion> >& fbc,
		       XMLReader& state_info_xml, 
		       const string& state_info_path,
		       const LinearOperator<LatticeFermion>& H);

    
    //! Create from OverlapStateInfo Structure
    const OverlapConnectState*
    createOverlapState(const multi1d<LatticeColorMatrix>& u_,
		       const FermBC< multi1d<LatticeFermion> >& fbc,
		       const OverlapStateInfo& state_info,
		       const LinearOperator<LatticeFermion>& H);
    

  }; // namespace OverlapConnectStateEnv

}; // namespace Chroma 

#endif
