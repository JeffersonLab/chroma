#ifndef FERM_HAMILTONIAN_H
#define FERM_HAMILTONIAN_H

#include "chromabase.h"
#include "fermact.h"
#include "update/molecdyn/abs_hamiltonian.h"
#include "actions/ferm/fermacts/prec_wilson_fermact_w.h"
#include "io/param_io.h"

template<typename GA, typename FA>
class TwoFlavorDegenFermHamiltonian : 
  public ExactAbsFermHamiltonian< multi1d<LatticeColorMatrix>,
				  multi1d<LatticeColorMatrix>,
				  LatticeFermion >
{
public: 

  // Constructor
  // A reference to the Gauge Action
  // A reference to an array of Fermion Actions
  // A reference to an array of Inverter Params
  TwoFlavorDegenFermHamiltonian(const GA& S_g_,
				const FA& S_f_,
				const InvertParam_t& inv_params_) : S_g(S_g_), S_f(S_f_), inv_params(inv_params_)  {}

  // Copy
  TwoFlavorDegenFermHamiltonian(const TwoFlavorDegenFermHamiltonian<GA,FA>& H) : S_g(H.S_g), S_f(H.S_f), inv_params(H.inv_params) {}

  // Clone
  virtual TwoFlavorDegenFermHamiltonian<GA,FA>* clone(void) const { 
    return new TwoFlavorDegenFermHamiltonian<GA,FA>(*this);
  }



  // Apply Gauge boundaries
  virtual void applyQBoundary(multi1d<LatticeColorMatrix>& q) const {
    S_g.getGaugeBC().modify(q);
  }

  virtual void applyPBoundary(multi1d<LatticeColorMatrix>& p) const {
    S_g.getGaugeBC().zero(p);
  }

  virtual void dsdq(const AbsPFFieldState<multi1d<LatticeColorMatrix>,
		    multi1d<LatticeColorMatrix>,LatticeFermion >& state,
		    multi1d<LatticeColorMatrix>&F) const {

    multi1d<LatticeColorMatrix> F_g(Nd);
    multi1d<LatticeColorMatrix> F_f(Nd);

    for(int mu=0; mu < Nd; mu++) { 
      F_g[mu] = zero;
      F_f[mu] = zero;
    }

    // Make a connect state with gauge boundaries
    const Handle< const ConnectState > g_bc_state(S_g.createState(state.getQ())); 
    // Get Gauge piece of force
    S_g.dsdu(F_g, g_bc_state);
    
    // Now I need to find X=(M^M)^{-1} \phi
    const LatticeFermion phi = state.getPhi()[0];
    LatticeFermion X = zero;
    
    // Make a connect state with fermionic boundaries
    const Handle< const ConnectState > f_bc_state(S_f.createState(state.getQ()));
    
    // Why a switch statement here? Later it may be useful for wilsoniums
    // to do a two stage BiCGStab ...
    switch(inv_params.invType) { 
    case CG_INVERTER:
      {
	// Get M
	Handle< const LinearOperator<LatticeFermion> > M(S_f.linOp(f_bc_state));
	
	// Get the subset for the linop
	const OrderedSubset& s = M->subset();

	int n_count=0;
	
	// Do the inversion M^{dag}M
	InvCG2(*M, phi, X, inv_params.RsdCG, inv_params.MaxCG, 
	       n_count);
	  
      }
      break;
    default:
      QDPIO::cerr << "Unsupported inverter type " << endl;
      QDP_abort(1);
      break;
    }

   
    // Now I can call S_f.dsdu
    S_f.dsdu(F_f, f_bc_state, X);

 

    // Add gauge and fermion contributions.
    for(int mu = 0; mu < Nd; mu++) { 
      F[mu] = F_g[mu] + F_f[mu];
    }
  }


  virtual Double mesKE(const AbsPFFieldState<multi1d<LatticeColorMatrix>,
		                             multi1d<LatticeColorMatrix>, 
                                             LatticeFermion>& state) const {
    // Extract momenta from the state
    // Once momenta have been refreshed, they are zeroed
    // on the gauge boundaries as requires so no need to mess
    // with gauge boundaries...
    const multi1d<LatticeColorMatrix>& mom = state.getP();

    // Square it up
    Double p_mom_sq = Double(0);


    // OK Here I am not sure about Sets/Subsets. GaugeAction has a 
    // GetSet() function but it is not a "subset" so do I just play 
    // on the whole lattice? I may need to insert some subsetting code
    // in here.
    // But for now just do the do...

    // This bit here is stolen from SZIN's MesMom()
    for(int mu = 0; mu < Nd; mu++) { 
      p_mom_sq += norm2(mom[mu]);
    }
    
    return p_mom_sq;
  }


  // Default -- Energy is: sum_phi phi^{dag} (MdagM)^{-1} phi
  // Note that this is particular to this Hamiltonian and NOT the action
  virtual Double mesFE(const AbsPFFieldState<multi1d<LatticeColorMatrix>,
		                             multi1d<LatticeColorMatrix>, 
		                             LatticeFermion>& state) const {
    
    // Get the Gauge Field
    const multi1d<LatticeColorMatrix>& u = state.getQ();
    Double ferm_energy = Double(0);
      
    // Loop over all the pseudofermion fields
    const LatticeFermion& phi = state.getPhi()[0];
    
    LatticeFermion psi = zero;
    
    // Make a connect state with fermionic boundaries
    const Handle< const ConnectState > g_state(S_f.createState(u));
      
    switch(inv_params.invType) { 
    case CG_INVERTER:
      {
	// Get M
	Handle< const LinearOperator<LatticeFermion> > M(S_f.linOp(g_state));
	
	// Get the subset for the linop
	const OrderedSubset& s = M->subset();

	int n_count=0;
	
	// Do the inversion M^{dag}M
	InvCG2(*M, phi, psi, inv_params.RsdCG, inv_params.MaxCG, 
	       n_count); 
	  
	// QDPIO::cout << "MesFE: n_count = " << n_count << endl;


	// Now get <phi, (MdagM)^{-1} phi> = <phi, psi>
 	ferm_energy = innerProductReal(phi,psi,s);
	
      }
      break;
    default:
      QDPIO::cerr << "Unsupported inverter type " << endl;
      QDP_abort(1);
      break;
    }
    return ferm_energy;
  }

  virtual Double mesGE(const AbsPFFieldState<multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix>, LatticeFermion>& state) const {

    Handle< const ConnectState > g_state(S_g.createState(state.getQ()));
    Double ret_val = S_g.S(g_state);
    return ret_val;
  }

 
  virtual ~TwoFlavorDegenFermHamiltonian(void) {};

  virtual const FA& getFermAct(const int i) const {
    return S_f;
  }

protected:

  
  GA S_g;
  FA S_f;
  InvertParam_t inv_params;
};




#endif
