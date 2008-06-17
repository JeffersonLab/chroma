// -*- C++ -*-
// $Id: eigen_state.h,v 1.2 2008-06-17 20:50:55 edwards Exp $
/*! \file
 *  \brief Eigenstate reader
 */

#ifndef eigen_state_h
#define eigen_state_h

#include "named_obj.h"
#include "state.h"
#include "chromabase.h"
#include "util/ferm/eigeninfo.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{

  //! Eigen-state holder
  /*! @ingroup fermstates */
  class EigenConnectState : public FermState<LatticeFermion, 
			                     multi1d<LatticeColorMatrix>, 
			                     multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Main constructor
    EigenConnectState(Handle< FermBC<T,P,Q> > fbc_, 
		      const multi1d<LatticeColorMatrix>& u_) : 
      fbc(fbc_), u(u_)
    {
      QDPIO::cout << "Calling EigenConnectState constructor with no IDs" << endl;
      fbc->modify(u);
      eigen_info_id = "";
      Neig = 0;
      dummy_evecs.resize(0);
      dummy_evals.resize(0);
    }

    EigenConnectState(Handle< FermBC<T,P,Q> > fbc_, 
		      const multi1d<LatticeColorMatrix>& u_,
		      std::string eigen_info_id_) :
      fbc(fbc_), u(u_)
    {
      fbc->modify(u);
      eigen_info_id = eigen_info_id_;
      QDPIO::cout << "Creating EigenConnectState using eigen_info_id :" << eigen_info_id << endl << flush ;

      Neig = TheNamedObjMap::Instance().getData< EigenInfo<T> >(eigen_info_id).getEvalues().size();
      dummy_evecs.resize(0);
      dummy_evals.resize(0);

    }

    ~EigenConnectState() {}

    const multi1d<LatticeColorMatrix>& getLinks() const { 
      return u;
    }

    //! Return the eigenvalues
    multi1d<Real>& getEvalues() {
      if (Neig == 0) { 
	return dummy_evals;
      }
      
      return TheNamedObjMap::Instance().getData< EigenInfo<T> >(eigen_info_id).getEvalues();
    }

    //! Return the eigenvalues
    const multi1d<Real>& getEvalues() const {
      if (Neig == 0) { 
	return dummy_evals;
      }
      
      return TheNamedObjMap::Instance().getData< EigenInfo<T> >(eigen_info_id).getEvalues();
    }

    multi1d<LatticeFermion>& getEvectors() { 
      if (Neig == 0) { 
	return dummy_evecs;
      }

      return TheNamedObjMap::Instance().getData< EigenInfo<T> >(eigen_info_id).getEvectors();
      
    }

    const multi1d<LatticeFermion>& getEvectors() const { 
      if (Neig == 0) { 
	return dummy_evecs;
      }

      return TheNamedObjMap::Instance().getData< EigenInfo<T> >(eigen_info_id).getEvectors();
      
    }

    Real& getLargest() {
      if (Neig == 0) { 
	QDPIO::cerr<< "Attempt to call getEvalues() on state with no e-values" << endl;
	QDP_abort(1);
      }
      
      return TheNamedObjMap::Instance().getData< EigenInfo<T> >(eigen_info_id).getLargest();      
    }

    const Real& getLargest() const {
      if (Neig == 0) { 
	QDPIO::cerr<< "Attempt to call getEvalues() on state with no e-values" << endl;
	QDP_abort(1);
      }
      
      return TheNamedObjMap::Instance().getData< EigenInfo<T> >(eigen_info_id).getLargest();      
    }
  
    int getNEig() const { 
      return Neig;
    }

    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

    //! Return the gauge BC object for this state
    /*! This is to help the optimized linops */
    Handle< FermBC<T,P,Q> > getFermBC() const {return fbc;}

  protected:
    //! Hide default constructor
    EigenConnectState() {}
    void operator=(const EigenConnectState&) {}

  private:
    Handle< FermBC<T,P,Q> > fbc;
    multi1d<LatticeColorMatrix> u;
    std::string eigen_info_id;
    int Neig;
    multi1d<LatticeFermion> dummy_evecs;
    multi1d<Real> dummy_evals;
  };



  //! Create a simple ferm connection state
  /*! @ingroup fermstates
   *
   * This is a factory class for producing a connection state
   */
  class CreateEigenConnectState : public CreateFermState<LatticeFermion, 
			                                 multi1d<LatticeColorMatrix>, 
			                                 multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Full constructor
    CreateEigenConnectState(Handle< FermBC<T,P,Q> > fbc_) : fbc(fbc_) {}

    //! Destructor
    ~CreateEigenConnectState() {}
   
    //! Construct a ConnectState
    EigenConnectState* operator()(const Q& q) const
      {
	return new EigenConnectState(fbc, q);
      }

    //! Return the ferm BC object for this state
    const FermBC<T,P,Q>& getBC() const {return *fbc;}

  private:
    CreateEigenConnectState() {}  // hide default constructur
    void operator=(const CreateEigenConnectState&) {} // hide =

  private:
    Handle< FermBC<T,P,Q> >  fbc;
  };



};



#endif
