// -*- C++ -*-

/*! @file 
 *  @brief Hyp field state for hyp links and a creator
 *
 *  Holds gauge fields at various smearing levels and caches some
 *  auxiliary things
 */

#ifndef _hyp_fermstate_h
#define _hyp_fermstate_h

#include "state.h"
#include "create_state.h"
#include "actions/ferm/fermstates/hyp_fermstate_params.h"
#include "util/gauge/hyp_utils.h"

namespace Chroma 
{
  /*! @ingroup fermstates */
  namespace CreateHypFermStateEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }


  //! Hyp field state
  /*! @ingroup fermstates
   *
   * Holds a hyp smeared state
   */
  template< typename T, typename P, typename Q>
  class HypFermState : public FermState<T,P,Q>
  {
  public: 
    // Typedefs to save typing
    //! Constructor only from a parameter structure
    // Constructor only from a parameter structure
    HypFermState(Handle< FermBC<T,P,Q> > fbc_, 
		   const HypFermStateParams& p_,
		   const Q& u_)
    {
      START_CODE();
      create(fbc_, p_, u_);
      END_CODE();
    }

    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

    //! Destructor is automagic
    virtual ~HypFermState() {}

    //! Return FAT Linke
    const Q& getLinks() const 
    {
      // This has BC-s applied already in the construction
      return fat_links_with_bc;
    }

    const Q& getThinLinks() const 
    {
      return smeared_links[0]; //return thin_links;
    }

    /* Recurse the thick force to compute the thin force. */
    /* I pull this out because SLIC like actions may want to do 
       this without slapping on the final list of thin U's 
       F_thin is resize internally */
    void fatForceToThin(const P& F_fat, P& F_thin) const
    {
      START_CODE();
      
      F_thin.resize(Nd);
      F_thin = F_fat;
      
      // Undo antiperiodic BC-s / force fixed BCs - this should really be unmodify
      // but essentially it works OK for everything except maybe twisted bc-s
      fbc->modify(F_thin);
      
      // Zero out fixed BCs
      fbc->zero(F_thin);
      
      // Now if the state is smeared recurse down.      
      for(int level=params.n_smear; level > 0; level--) {
	
	Hyping::deriv_recurse(F_thin, params.smear_in_this_dirP, 
                              params.BlkMax,
                              params.BlkAccu,
                              smeared_links[level-1]);
	
	fbc->zero(F_thin);
	
      }
      
      END_CODE();
    }
    

    /* This recurses the force and slaps on the gauge piece at the end.
       I am coming around to realising that the gauge stuff ought not 
       be slapped on here but maybe somewhere in the MC. Consider a SLIC
       force. There the operator would appply the fatForceToThin to 
       the fat force and just do the normal deriv for the thin liks
       and this thing can be pulled out into the dsdq methods.
       For now, the default behaviour is to recurse the force down 
       here but I am planning ahead to a refactor. */
    virtual void deriv(P& F) const 
    {
      START_CODE();
      
      P F_tmp(Nd);

      // Function resizes F_tmp
      fatForceToThin(F,F_tmp);
      
      // Multiply in by the final U term to close off the links
      for(int mu=0; mu < Nd; mu++) { 
	F[mu] = (smeared_links[0])[mu]*F_tmp[mu];
      }

      END_CODE();
    }


  private:

    // Hide default constructor
    HypFermState(){}


    // create function
    void create(Handle< FermBC<T,P,Q> > fbc_,
		const HypFermStateParams& p_,
		const Q& u_)
    { 
      START_CODE();
      
      fbc    =  fbc_ ;
      params =  p_   ;

      // Allocate smeared and thin links
      smeared_links.resize(params.n_smear + 1);
      for(int i=0; i <= params.n_smear; i++) smeared_links[i].resize(Nd);

      // Copy thin links into smeared_links[0]
      for(int mu=0; mu < Nd; mu++) (smeared_links[0])[mu] = u_[mu];

      // Apply BC
      if( fbc->nontrivialP() ) fbc->modify( smeared_links[0] );    

      // Iterate up the smearings
      for(int i=1; i <= params.n_smear; i++) {
	
	Hyping::smear_links(smeared_links[i-1], smeared_links[i], 
                            params.smear_in_this_dirP, 
                            params.alpha1,
                            params.alpha2,
                            params.alpha3,
                            params.BlkMax,
                            params.BlkAccu);
        
	if( fbc->nontrivialP() ) {
	  fbc->modify( smeared_links[i] );    
	}	
      }

      // ANTIPERIODIC BCs only -- modify only top level smeared thing
      fat_links_with_bc.resize(Nd);
      fat_links_with_bc = smeared_links[params.n_smear];
      fbc->modify(fat_links_with_bc);

      END_CODE();
    }
    

  private:
    Handle< FermBC<T,P,Q> >  fbc;

    // smeared_links[0] are the thin links smeared_links[params.n_smear] 
    // are the smeared links.
    multi1d< Q > smeared_links;
    Q fat_links_with_bc;
    
    HypFermStateParams  params;

  }; // End class



  //! Create a hyp ferm connection state
  /*! @ingroup fermstates
   *
   * This is a factory class for producing a connection state
   */
  template<typename T, typename P, typename Q >
  class CreateHypFermState : public CreateFermState<T, P, Q>
  {
  public: 
    // Typedefs to save typing
    
    //! Full constructor
    CreateHypFermState(Handle< FermBC<T,P,Q> > fbc_,
			 const HypFermStateParams& p_) : 
      fbc(fbc_), params(p_) {}
    
    //! Destructor
    ~CreateHypFermState() {}
    
    //! Construct a ConnectState
    HypFermState<T,P,Q>* operator()(const Q& q) const
    {
      return new HypFermState<T,P,Q>(fbc, params, q);
    }
    
    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

  private:
    CreateHypFermState() {}  // hide default constructur
    void operator=(const CreateHypFermState&) {} // hide =

  private:
    Handle< FermBC<T,P,Q> >  fbc;
    HypFermStateParams  params;
  };


}
#endif
