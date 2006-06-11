#include "chromabase.h"

#include "actions/ferm/qprop/asqtad_cps_wrapper_qprop.h"

namespace Chroma 
{

  //! Constructor
  /*!
  // Keeping the same interface as for the ordinary staggered 
  // qprop...
  //
  // But the M_ and A_ linop handles are no longer kept
  // (are ignored) -- is there a nice way around this ? 
  // Perhaps not
  */
  AsqtadCPSWrapperQprop::AsqtadCPSWrapperQprop(EvenOddStaggeredTypeFermAct<LatticeStaggeredFermion, 
					       multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& S_,
					       Handle< FermState<LatticeStaggeredFermion,
					       multi1d<LatticeColorMatrix>,
					       multi1d<LatticeColorMatrix> > state_, 
					       const InvertParam_t& invParam_) :
    Mass(S_.getQuarkMass()), invParam(invParam_), state(state_)
  {
    
    // Do nothing for now. Zbysh fill in as needed
    
    
  }
  
  
  //! Destructor may not well be automatic
  AsqtadCPSWrapperQprop::~AsqtadCPSWrapperQprop() 
  {
    
    // Zbysh -- do the destruction here.
  }


  SystemSolverResults_t
  AsqtadCPSWrapperQprop::operator() (LatticeStaggeredFermion& psi, 
				     const LatticeStaggeredFermion& chi) const
  {
    SystemSolverResults_t res;

    // Here is how to get a the mass
    QDPIO::cout << "Mass is " << Mass << endl;

    // Here is how to get at the gauge links:
    const multi1d<LatticeColorMatrix>& u = state->getLinks();
    
    
    // Do what thou wouldst here
    QDPIO::cout << "in Asqtad_CPS_qprop_wrapper::operator() stub" << endl;
    QDPIO::cout << "Rewrite this bit please" << endl;

    res.n_count = 0;
    psi = chi;
    
    return res;
  }

}
