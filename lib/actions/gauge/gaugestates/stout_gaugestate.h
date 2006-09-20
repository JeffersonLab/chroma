// -*- C++ -*-
// $Id: stout_gaugestate.h,v 1.2 2006-09-20 20:28:01 edwards Exp $

/*! @file
 * @brief Stout gauge state and a creator
 */

#ifndef __stout_gaugestate_h__
#define __stout_gaugestate_h__

#include "state.h"
#include "create_state.h"
#include "handle.h"

#include "actions/ferm/fermstates/stout_fermstate_params.h"
#include "util/gauge/stout_utils.h"

namespace Chroma
{

  /*! @ingroup gaugestates */
  namespace CreateStoutGaugeStateEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }


  //! Stout version of GaugeState 
  /*! @ingroup gaugestates
   *
   * Only needs to hold a gauge field and gauge bc
   */
  template<typename P, typename Q>
  class StoutGaugeState : public GaugeState<P,Q>
  {
  public:
    //! Full constructor
    StoutGaugeState(Handle< GaugeBC<P,Q> > gbc_,
		    const StoutFermStateParams& p_,
		    const Q& q_)
      {
	START_CODE();
	create(gbc_, p_, q_);
	END_CODE();
      }
    
    //! Return the gauge BC object for this state
   const GaugeBC<P,Q>& getBC() const {return *gbc;}
    
    //! Destructor
    virtual ~StoutGaugeState() {}

    //! Return the link fields needed in constructing linear operators
    const Q& getLinks() const {
      return smeared_links[ params.n_smear ];
    }

    const Q& getThinLinks() const {
      return smeared_links[0];
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
      
      // Zero out fixed BCs
      // gbc->zero(F_thin);
      
      // Now if the state is smeared recurse down.
      
      for(int level=params.n_smear; level > 0; level--) {
	QDPIO::cout << "Recursing level " << level << " to level " << level -1 << endl;
	Stouting::deriv_recurse(F_thin, params.smear_in_this_dirP, params.rho, smeared_links[level-1]);
	
	// gbc->zero(F_thin);
	
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
    void deriv(P& F) const 
    {
      START_CODE();

      QDPIO::cout << "WARNING: StoutGaugeState::deriv() is never called in current scheme of things -- the gauge monomial doesn't know about it" << endl;
      
      P F_tmp(Nd);
      
      // Function resizes F_tmp
      fatForceToThin(F,F_tmp);
      
      
      // Multiply in by the final U term to close off the links
      for(int mu=0; mu < Nd; mu++) { 
	F[mu] = (smeared_links[0])[mu]*F_tmp[mu];
      }
      
      END_CODE();
    }
    

    void create(Handle< GaugeBC<P,Q> > gbc_,
		const StoutFermStateParams& p_,
		const Q& u_)
    { 
      START_CODE();
      
      gbc = gbc_;
      params = p_;
    
      // Allocate smeared and thin links
      smeared_links.resize(params.n_smear + 1);
      for(int i=0; i <= params.n_smear; i++) { 
	smeared_links[i].resize(Nd);
      }
      
      
      // Copy thin links into smeared_links[0]
      for(int mu=0; mu < Nd; mu++) { 
	(smeared_links[0])[mu] = u_[mu];
      }
      
      // This is different from fermions.
      // I apply the BC-s to the gauge fields
      // at the bottom. 
      gbc->modify( smeared_links[0] );    
      
      // Iterate up the smearings
      for(int i=1; i <= params.n_smear; i++) {
	Stouting::smear_links(smeared_links[i-1], smeared_links[i], params.smear_in_this_dirP, params.rho);

	// If the gauge BC's are nontrivial 
	// apply at every level - eg SF BC's
	if( gbc->nontrivialP() ) {
	  gbc->modify( smeared_links[i] );    
        }
	
      }
      
      END_CODE();
    }

  private:
    StoutGaugeState() {}  // hide default constructur
    void operator=(const StoutGaugeState&) {} // hide =

  private:
    Handle< GaugeBC<P,Q> >  gbc;
    multi1d<Q> smeared_links;
    StoutFermStateParams params;
  };



  //! Create a stout gauge connection state
  /*! @ingroup gaugestates
   *
   * This is a factory class for producing a connection state
   */
  template<typename P, typename Q>
  class CreateStoutGaugeState : public CreateGaugeState<P,Q>
  {
  public:
    //! Full constructor
    CreateStoutGaugeState(Handle< GaugeBC<P,Q> > gbc_,
			  const StoutFermStateParams& p_) : gbc(gbc_), params(p_) {}

    //! Destructor
    ~CreateStoutGaugeState() {}
   
    //! Construct a ConnectState
    StoutGaugeState<P,Q>* operator()(const Q& q) const
      {
	return new StoutGaugeState<P,Q>(gbc, params, q);
      }

    //! Return the gauge BC object for this state
    const GaugeBC<P,Q>& getBC() const {return *gbc;}

  private:
    CreateStoutGaugeState() {}  // hide default constructur
    void operator=(const CreateStoutGaugeState&) {} // hide =

  private:
    Handle< GaugeBC<P,Q> >  gbc;
    StoutFermStateParams params;
  };



}


#endif
