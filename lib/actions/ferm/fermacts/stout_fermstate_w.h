// -*- C++ -*-
// $Id: stout_fermstate_w.h,v 1.6 2006-08-15 13:17:24 bjoo Exp $

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
    ~StoutFermState() {}

    //! Return FAT Linke
    const multi1d<LatticeColorMatrix>& getLinks() const 
    {
      // This has BC-s applied already in the construction
      return smeared_links[params.n_smear];
    }

    const multi1d<LatticeColorMatrix>& getThinLinks() const 
    {
      // This has the BC-s applied already
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
    void deriv(multi1d<LatticeColorMatrix>& F) const {

      multi1d<LatticeColorMatrix> F_tmp;

      // Function resizes F_tmp
      fatForceToThin(F,F_tmp);
      
      // Multiply in by the final U term to close off the links
      for(int mu=0; mu < Nd; mu++) { 
	F[mu] = (smeared_links[0])[mu]*F_tmp[mu];
      }
    }

    //! given field U, construct the staples into C, form Q and Q^2 and compute
    //  c0 and c1
    void getQsandCs(const multi1d<LatticeColorMatrix>& u, 
		    LatticeColorMatrix& Q, 
		    LatticeColorMatrix& QQ,
		    LatticeColorMatrix& C, 
		    LatticeDouble& c0,
		    LatticeDouble& c1,
		    int mu) const;

    //! Given c0 and c1 compute the f-s and b-s
    //! Only compute b-s if do_bs is set to true (default)
    void getFsAndBs(const LatticeDouble& c0,
		    const LatticeDouble& c1,
		    multi1d<LatticeDComplex>& f,
		    multi1d<LatticeDComplex>& b1,
		    multi1d<LatticeDComplex>& b2,
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

} // end namespace

#endif
