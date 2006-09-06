// -*- C++ -*-
// $Id: stout_fermstate_w.h,v 1.10 2006-09-06 13:18:22 bjoo Exp $

/*! @file 
 *  @brief Stout field state for stout links and a creator
 *
 *  Holds gauge fields at various smearing levels and caches some
 *  auxiliary things
 */

#ifndef _stout_fermstate_h
#define _stout_fermstate_h

#include "state.h"
#include "create_state.h"
#include "actions/ferm/fermacts/stout_fermstate_params.h"

namespace Chroma 
{
  /*! @ingroup fermacts */
  namespace CreateStoutFermStateEnv 
  { 
    extern const std::string name;
    extern const bool registered;
  }
  
  namespace StoutLinkTimings { 
    double getForceTime();
    double getSmearingTime();
    double getFunctionsTime();
  };

  //! Stout field state
  /*! @ingroup fermacts
   *
   * Holds a stout smeared state
   */
  class StoutFermState : public FermState<LatticeFermion,
			 multi1d<LatticeColorMatrix>,
			 multi1d<LatticeColorMatrix> >
  {
  public: 
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Constructor only from a parameter structure
    StoutFermState(Handle< FermBC<T,P,Q> > fbc_, 
		   const StoutFermStateParams& p_,
		   const multi1d<LatticeColorMatrix>& u_);


    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

    //! Destructor is automagic
    virtual ~StoutFermState() {}

    //! Return FAT Linke
    const multi1d<LatticeColorMatrix>& getLinks() const 
    {
      // This has BC-s applied already in the construction
      return fat_links_with_bc;
    }

    const multi1d<LatticeColorMatrix>& getThinLinks() const 
    {
      return smeared_links[0];
    }

    /* Recurse the thick force to compute the thin force. */
    /* I pull this out because SLIC like actions may want to do 
       this without slapping on the final list of thin U's 
       F_thin is resize internally */
    void fatForceToThin(const multi1d<LatticeColorMatrix>& F_fat, multi1d<LatticeColorMatrix>& F_thin) const;

    /* This recurses the force and slaps on the gauge piece at the end.
       I am coming around to realising that the gauge stuff ought not 
       be slapped on here but maybe somewhere in the MC. Consider a SLIC
       force. There the operator would appply the fatForceToThin to 
       the fat force and just do the normal deriv for the thin liks
       and this thing can be pulled out into the dsdq methods.
       For now, the default behaviour is to recurse the force down 
       here but I am planning ahead to a refactor. */
    virtual void deriv(multi1d<LatticeColorMatrix>& F) const;

   
    //! given field U, construct the staples into C, form Q and Q^2 and compute
    //  c0 and c1
    void getQsandCs(const multi1d<LatticeColorMatrix>& u, 
		    LatticeColorMatrix& Q, 
		    LatticeColorMatrix& QQ,
		    LatticeColorMatrix& C, 
		    int mu) const;

    //! Given c0 and c1 compute the f-s and b-s
    //! Only compute b-s if do_bs is set to true (default)
    void getFsAndBs(const LatticeColorMatrix& Q,
		    const LatticeColorMatrix& QQ,
		    multi1d<LatticeComplex>& f,
		    multi1d<LatticeComplex>& b1,
		    multi1d<LatticeComplex>& b2,
		    bool do_bs=true) const;

  private:

    // Hide default constructor
    StoutFermState(){}

    // Do the smearing from level i to level i+1
    void smear_links(const multi1d<LatticeColorMatrix>& current,
		     multi1d<LatticeColorMatrix>& next);


    // Do the force recursion from level i+1, to level i
    void deriv_recurse(multi1d<LatticeColorMatrix>&  F,
		       const int level) const;


    // create function
    void create(Handle< FermBC<T,P,Q> > fbc_,
		const StoutFermStateParams& p_,
		const multi1d<LatticeColorMatrix>& u_);


  private:
    Handle< FermBC<T,P,Q> >  fbc;

    // smeared_links[0] are the thin links smeared_links[params.n_smear] 
    // are the smeared links.
    multi1d< multi1d<LatticeColorMatrix> > smeared_links;
    multi1d<LatticeColorMatrix> fat_links_with_bc;


    StoutFermStateParams  params;
  }; // End class


  //! Create a stout ferm connection state
  /*! @ingroup fermacts
   *
   * This is a factory class for producing a connection state
   */
  class CreateStoutFermState : public CreateFermState<LatticeFermion,
			       multi1d<LatticeColorMatrix>,
			       multi1d<LatticeColorMatrix> >
  {
  public: 
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;
    
    //! Full constructor
    CreateStoutFermState(Handle< FermBC<T,P,Q> > fbc_,
			 const StoutFermStateParams& p_) : 
      fbc(fbc_), params(p_) {}
    
    //! Destructor
    ~CreateStoutFermState() {}
    
    //! Construct a ConnectState
    StoutFermState* operator()(const Q& q) const
    {
      return new StoutFermState(fbc, params, q);
    }
    
    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

  private:
    CreateStoutFermState() {}  // hide default constructur
    void operator=(const CreateStoutFermState&) {} // hide =

  private:
    Handle< FermBC<T,P,Q> >  fbc;
    StoutFermStateParams  params;
  };


  /* SLIC Fermions */

  /*! @ingroup fermacts */
  namespace CreateSLICFermStateEnv 
  { 
    extern const std::string name;
    extern const bool registered;
  }

  //! SLIC (Stout Link Irrelevant Clover ferm connection state
  /*! @ingroup fermacts
   *
   * This ferm state is for use in SLIC Fermions ONLY (or their ilk)
   * It inherits directly from stout ferm state and overrides the 
   * deriv function. In the SLIC case it is the duty of the linop 
   * to thin the fat links, so here the deriv() only multiplies by 
   * the thin links. CAVEAT: The way stout links currently work
   * The fermion force for the stout part must remove the fermionic
   * boundaries - as the thin links don't have these applied. THe 
   * thin links do have gauge boundaries (a la schroedinger) applied 
   * 
   *  - NB This state inherits data from its super. This is not necessarily 
   *    a good thing. However we don't have private data of our own that is 
   *    not in the superclass and I only access this read only in deriv()
   *    using a public function.
   */
  class SLICFermState : public StoutFermState {
  public:
    // Initialise the super
    SLICFermState(Handle< FermBC<T,P,Q> > fbc_, 
		  const StoutFermStateParams& p_,
		  const multi1d<LatticeColorMatrix>& u_) : StoutFermState(fbc_,p_, u_) {}

    ~SLICFermState() {}
   
    void deriv(multi1d<LatticeColorMatrix>& F) const {

      // Multiply in by the final U term to close off the links
      for(int mu=0; mu < Nd; mu++) { 
	F[mu] = getThinLinks()[mu]*F[mu];
      }
    }
  };

  //! Create a SLIC ferm connection state
  /*! @ingroup fermacts
   *
   * This is a factory class for producing a connection state
   */
  class CreateSLICFermState : public CreateFermState<LatticeFermion,
			      multi1d<LatticeColorMatrix>,
			      multi1d<LatticeColorMatrix> >
  {
  public: 
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    CreateSLICFermState(Handle< FermBC<T,P,Q> > fbc_,
			const StoutFermStateParams& p_) : 
      fbc(fbc_), params(p_) {}

    //! Destructor
    ~CreateSLICFermState() {}
   
    //! Construct a ConnectState
    StoutFermState* operator()(const Q& q) const
      {
	return new SLICFermState(fbc, params, q);
      }

    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

  private:
    CreateSLICFermState() {}  // hide default constructur
    void operator=(const CreateSLICFermState&) {} // hide =

  private:
    Handle< FermBC<T,P,Q> >  fbc;
    StoutFermStateParams  params;
  };

} // end namespace

#endif
