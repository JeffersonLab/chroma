// -*- C++ -*-
// $Id: stout_fermstate_w.h,v 1.5 2007-11-28 21:11:08 bjoo Exp $

/*! @file 
 *  @brief Hex field state for stout links and a creator
 *
 *  Holds gauge fields at various smearing levels and caches some
 *  auxiliary things
 */

#ifndef _hex_fermstate_h
#define _hex_fermstate_h

#include "state.h"
#include "create_state.h"
#include "actions/ferm/fermstates/hex_fermstate_params.h"
#include "util/gauge/stout_utils.h"
#include "meas/smear/hex_smear.h"

namespace Chroma 
{
  /*! @ingroup fermstates */
  namespace CreateHexFermStateEnv 
  { 
    extern const std::string name;
    bool registerAll();
  }


  //! Hex field state
  /*! @ingroup fermstates
   *
   * Holds a hex smeared state
   */

  template< typename T, typename P, typename Q>
  class HexFermState : public FermState<T,P,Q>
  {
  public: 
    // Typedefs to save typing
    //! Constructor only from a parameter structure
  // Constructor only from a parameter structure
    HexFermState(Handle< FermBC<T,P,Q> > fbc_, 
		   const HexFermStateParams& p_,
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
    virtual ~HexFermState() {}

    //! Return FAT Linke
    const Q& getLinks() const 
    {
      // This has BC-s applied already in the construction
      return fat_links_with_bc;
    }

    const Q& getThinLinks() const 
    {
      return thin_links;
    }

    /* Only an inverter at the moment */
    void fatForceToThin(const P& F_fat, P& F_thin) const
    {
      START_CODE();
      
      
      END_CODE();
    }
    

    /* Not implemneted   **/

    virtual void deriv(P& F) const 
    {
      START_CODE();
      
      END_CODE();
    }


  private:

    // Hide default constructor
    HexFermState(){}


    // create function
    void create(Handle< FermBC<T,P,Q> > fbc_,
		const HexFermStateParams& p_,
		const Q& u_)
    { 
      START_CODE();
      
      fbc    =  fbc_ ;
      params =  p_   ;
      
      // allocate some space
      thin_links.resize(Nd);
      smeared_links.resize(Nd);

      thin_links = u_; 

      Hex_Smear(u_, smeared_links, params.n_smear) ;

      
      if( fbc->nontrivialP() ) {
	fbc->modify( smeared_links );    
      }
      
      // ANTIPERIODIC BCs only -- modify only top level smeared thing
      fat_links_with_bc.resize(Nd)      ;
      fat_links_with_bc = smeared_links ;
      fbc->modify(fat_links_with_bc)    ;
      
      
      END_CODE();
    }
    

  private:
    Handle< FermBC<T,P,Q> >  fbc;
    
    Q  thin_links      ;
    Q  smeared_links      ;
    Q  fat_links_with_bc  ;
    
    HexFermStateParams  params;

  }; // End class



  //! Create a hex ferm connection state
  /*! @ingroup fermstates
   *
   * This is a factory class for producing a connection state
   */
  template<typename T, typename P, typename Q >
  class CreateHexFermState : public CreateFermState<T, P, Q>
  {
  public: 
    // Typedefs to save typing
    
    //! Full constructor
    CreateHexFermState(Handle< FermBC<T,P,Q> > fbc_,
			 const HexFermStateParams& p_) : 
      fbc(fbc_), params(p_) {}
    
    //! Destructor
    ~CreateHexFermState() {}
    
    //! Construct a ConnectState
    HexFermState<T,P,Q>* operator()(const Q& q) const
    {
      return new HexFermState<T,P,Q>(fbc, params, q);
    }
    
    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

  private:
    CreateHexFermState() {}  // hide default constructur
    void operator=(const CreateHexFermState&) {} // hide =

  private:
    Handle< FermBC<T,P,Q> >  fbc;
    HexFermStateParams  params;
  };


}
#endif
