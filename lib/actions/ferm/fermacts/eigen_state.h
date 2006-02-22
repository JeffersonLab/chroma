#ifndef eigen_state_h
#define eigen_state_h

#include "named_obj.h"
#include "state.h"
#include "chromabase.h"
#include "util/ferm/eigeninfo.h"
#include "meas/inline/io/named_objmap.h"
namespace Chroma {



  class EigenConnectState : public ConnectState {
  public:
    EigenConnectState(const multi1d<LatticeColorMatrix>& u_) 
    {
      QDPIO::cout << "Calling EigenConnectState constructor with no IDs" << endl;
      u = u_;
      eigen_info_id = "";
      Neig = 0;
      dummy_evecs.resize(0);
      dummy_evals.resize(0);
    }

    EigenConnectState(const multi1d<LatticeColorMatrix>& u_,
		      std::string eigen_info_id_) {
      u = u_;
      eigen_info_id = eigen_info_id_;
      QDPIO::cout << "Creating EigenConnectState using eigen_info_id :" << eigen_info_id << endl << flush ;

      Neig = TheNamedObjMap::Instance().getData<EigenInfo>(eigen_info_id).getEvalues().size();
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
      
      return TheNamedObjMap::Instance().getData<EigenInfo>(eigen_info_id).getEvalues();
      
    }

    //! Return the eigenvalues
    const multi1d<Real>& getEvalues() const {
      if (Neig == 0) { 
	return dummy_evals;
      }
      
      return TheNamedObjMap::Instance().getData<EigenInfo>(eigen_info_id).getEvalues();
      
    }

    multi1d<LatticeFermion>& getEvectors() { 
      if (Neig == 0) { 
	return dummy_evecs;
      }

      return TheNamedObjMap::Instance().getData<EigenInfo>(eigen_info_id).getEvectors();
      
    }

    const multi1d<LatticeFermion>& getEvectors() const { 
      if (Neig == 0) { 
	return dummy_evecs;
      }

      return TheNamedObjMap::Instance().getData<EigenInfo>(eigen_info_id).getEvectors();
      
    }

    Real& getLargest() {
      if (Neig == 0) { 
	QDPIO::cerr<< "Attempt to call getEvalues() on state with no e-values" << endl;
	QDP_abort(1);
      }
      
      return TheNamedObjMap::Instance().getData<EigenInfo>(eigen_info_id).getLargest();      
    }

    const Real& getLargest() const {
      if (Neig == 0) { 
	QDPIO::cerr<< "Attempt to call getEvalues() on state with no e-values" << endl;
	QDP_abort(1);
      }
      
      return TheNamedObjMap::Instance().getData<EigenInfo>(eigen_info_id).getLargest();      
    }
  
    int getNEig() const { 
      return Neig;
    }

  private:
    multi1d<LatticeColorMatrix> u;
    std::string eigen_info_id;
    int Neig;
    multi1d<LatticeFermion> dummy_evecs;
    multi1d<Real> dummy_evals;
  };




};



#endif
