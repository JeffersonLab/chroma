// -*- C++ -*-
// $Id: overlap_state.h,v 1.1 2006-09-19 17:53:37 edwards Exp $
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
  
  //! Overlap connection state
  /*! \ingroup fermstates */
  class OverlapConnectState : public FermState<LatticeFermion,
			      multi1d<LatticeColorMatrix>,
			      multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    typedef Real WordBase_t;
  
    //! Constructor with no eigenvalues
    OverlapConnectState(Handle< FermBC<T,P,Q> > fbc_,
			const multi1d<LatticeColorMatrix>& u_,  // gauge field
			const WordBase_t& approxMin_,          // epsilon
			const WordBase_t& approxMax_);         // approx max

    //! Constructor with e-values and e-vectors
    OverlapConnectState(Handle< FermBC<T,P,Q> > fbc_,
			const multi1d<LatticeColorMatrix>& u_,
			const multi1d<WordBase_t>& val_, 
			const multi1d<LatticeFermion>& vec_,
			const WordBase_t& val_max_,
			const WordBase_t& approxMin_,
			const WordBase_t& approxMax_);

    //---------------------------------------------
    // These versions of constructors deal with funky input. RGE would like
    // to get rid of them.
    
    //! Create a ConnectState with just the gauge fields
    OverlapConnectState(Handle< FermBC<T,P,Q> > fbc_,
			const multi1d<LatticeColorMatrix>& u_);
    
    //! Create a ConnectState with just the gauge fields, and a lower
    //!  approximation bound
    OverlapConnectState(Handle< FermBC<T,P,Q> > fbc_,
			const multi1d<LatticeColorMatrix>& u_,
			const Real& approxMin_);

    //! Create OverlapConnectState with eigenvalues/vectors
    OverlapConnectState(Handle< FermBC<T,P,Q> > fbc_,
			const multi1d<LatticeColorMatrix>& u_,
			const multi1d<Real>& lambda_lo_, 
			const multi1d<LatticeFermion>& evecs_lo_,
			const Real& lambda_hi_);

    //! Create a ConnectState out of XML
    OverlapConnectState(Handle< FermBC<T,P,Q> > fbc_,
			const multi1d<LatticeColorMatrix> u_,
			XMLReader& state_info_xml, 
			const string& state_info_path,
			const LinearOperator<LatticeFermion>& H);

    //! Create from OverlapStateInfo Structure
    OverlapConnectState(Handle< FermBC<T,P,Q> > fbc_,
			const multi1d<LatticeColorMatrix>& u_,
			const OverlapStateInfo& state_info,
			const LinearOperator<LatticeFermion>& H);
    //---------------------------------------------

    //! Copy constructor
    OverlapConnectState(const OverlapConnectState& a) : fbc(a.fbc), u(a.u), 
							eigVal(a.eigVal), eigVec(a.eigVec), 
							eigValMax(a.eigValMax), 
							approxMin(a.approxMin), 
							approxMax(a.approxMax) {}

    ~OverlapConnectState() {};

    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

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

    //! for now inherit the deriv operation

  protected:
    //! Constructor with no eigenvalues
    void init(Handle< FermBC<T,P,Q> > fbc_,
	      const multi1d<LatticeColorMatrix>& u_,  // gauge field
	      const WordBase_t& approxMin_,          // epsilon
	      const WordBase_t& approxMax_);         // approx max

    //! Constructor with e-values and e-vectors
    void init(Handle< FermBC<T,P,Q> > fbc_,
	      const multi1d<LatticeColorMatrix>& u_,
	      const multi1d<WordBase_t>& val_, 
	      const multi1d<LatticeFermion>& vec_,
	      const WordBase_t& val_max_,
	      const WordBase_t& approxMin_,
	      const WordBase_t& approxMax_);

    //! Create from OverlapStateInfo Structure
    void init(Handle< FermBC<T,P,Q> > fbc_,
	      const multi1d<LatticeColorMatrix>& u_,
	      const OverlapStateInfo& state_info,
	      const LinearOperator<LatticeFermion>& H);

  private:
    OverlapConnectState() {}  // hide default constructur
    void operator=(const OverlapConnectState&) {} // hide =

  private:
    Handle< FermBC<T,P,Q> > fbc;
    multi1d<LatticeColorMatrix> u;
    multi1d<WordBase_t>  eigVal;
    multi1d<LatticeFermion> eigVec;
    WordBase_t           eigValMax;
    WordBase_t           approxMin; 
    WordBase_t           approxMax;
  };



#if 0
  // NOTE: there should be a a creator here with some of the functionality
  // of the overlapstate moved inside. This might simply a tad bit the
  // creation. E.g., all the funky ways to create a state would be moved
  // here and the state itself only hold the basic stuff.
  class CreateOverlapConnectState : public CreateFermState<LatticeFermion,
				    multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> >
  {
  public:
    //! Constructor with no eigenvalues
    CreateOverlapConnectState();
  };
#endif

   
} // namespace Chroma 

#endif
