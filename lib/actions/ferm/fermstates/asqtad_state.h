// -*- C++ -*-
// $Id: asqtad_state.h,v 1.1 2006-09-19 17:53:36 edwards Exp $
/*! \file
 *  \brief Asqtad state
 */

#ifndef __asqtad_state_h__
#define __asqtad_state_h__


#include "handle.h"
#include "state.h"
#include "actions/ferm/linop/improvement_terms_s.h"


namespace Chroma 
{ 
  //! Basic "Connect State" for ASQTAD
  /*! 
   * \ingroup fermstates
   */
  class AsqtadConnectStateBase : public FermState<LatticeStaggeredFermion, 
				 multi1d<LatticeColorMatrix>, 
				 multi1d<LatticeColorMatrix> >
  {
  public: 

    virtual const multi1d<LatticeColorMatrix>& getFatLinks() const = 0;
    virtual const multi1d<LatticeColorMatrix>& getTripleLinks() const = 0;

    //! Return the gauge BC object for this state
    virtual const FermBC<LatticeStaggeredFermion, 
			 multi1d<LatticeColorMatrix>, 
			 multi1d<LatticeColorMatrix> >& getBC() const = 0;

    //! Return the ferm BC object for this state
    virtual Handle< FermBC<LatticeStaggeredFermion, 
			   multi1d<LatticeColorMatrix>, 
			   multi1d<LatticeColorMatrix> > > getFermBC() const = 0;

  };


  //! The actual Asqtad thing
  /*! 
   * \ingroup fermstates
   */
  class AsqtadConnectState : public AsqtadConnectStateBase
  {
  public:
    // Typedefs to save typing
    typedef LatticeStaggeredFermion      T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    typedef Real WordBase_t;
  
    //! Full Constructor
    // Make deep copies here
    AsqtadConnectState(Handle< FermBC<T,P,Q> > fbc_,
		       const multi1d<LatticeColorMatrix>& u_,
		       const multi1d<LatticeColorMatrix>& u_fat_,
		       const multi1d<LatticeColorMatrix>& u_triple_)
      : fbc(fbc_), u(u_), u_fat(u_fat_), u_triple(u_triple_)  { };

    ~AsqtadConnectState() {};

    //! Return the link fields needed in constructing linear operators
    const multi1d<LatticeColorMatrix>& getLinks() const { return u; }
    const multi1d<LatticeColorMatrix>& getFatLinks() const { return u_fat; }
    const multi1d<LatticeColorMatrix>& getTripleLinks() const { return u_triple; }

    //! Return the gauge BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}
   
    //! Return the ferm BC object for this state
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

  private:
    AsqtadConnectState() {}  // hide default constructur
    void operator=(const AsqtadConnectState&) {} // hide =

  private:
    Handle< FermBC<T,P,Q> > fbc;
    multi1d<LatticeColorMatrix> u;
    multi1d<LatticeColorMatrix> u_fat;
    multi1d<LatticeColorMatrix> u_triple;
  };


} // End Namespace Chroma


#endif
