// $Id: unprec_wilson_fermact_w.cc,v 1.22 2004-12-07 17:11:50 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/fermfactory_w.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace UnprecWilsonFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion>* createFermAct(Handle< FermBC<LatticeFermion> > fbc,
						     XMLReader& xml_in,
						     const std::string& path)
    {
      return new UnprecWilsonFermAct(fbc, UnprecWilsonFermActParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_WILSON";

    //! Register the Wilson fermact
    const bool registered = TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct);
  }


  //! Read parameters
  UnprecWilsonFermActParams::UnprecWilsonFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    if (paramtop.count("Mass") != 0) 
    {
      read(paramtop, "Mass", Mass);
      if (paramtop.count("Kappa") != 0) 
      {
	QDPIO::cerr << "Error: found both a Kappa and a Mass tag" << endl;
	QDP_abort(1);
      }
    }
    else if (paramtop.count("Kappa") != 0)
    {
      Real Kappa;
      read(paramtop, "Kappa", Kappa);
      Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
    }
    else
    {
      QDPIO::cerr << "Error: neither Mass or Kappa found" << endl;
      QDP_abort(1);
    }

    // There is always an aniso Param for wilson, so set it to default
    initHeader(anisoParam);

    //  Read optional aniso
    if (paramtop.count("AnisoParam") != 0) 
      read(paramtop, "AnisoParam", anisoParam);
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecWilsonFermActParams& param)
  {
    UnprecWilsonFermActParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator<LatticeFermion>*
  UnprecWilsonFermAct::linOp(Handle<const ConnectState> state) const
  {
    if (param.anisoParam.anisoP)
    {
      QDPIO::cerr << "UnprecWilsonFermAct::linOp - no aniso support" << endl;
      QDP_abort(1);
    }

    return new UnprecWilsonLinOp(state->getLinks(), param.Mass); 
  }

  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state    gauge field     	       (Read)
   */
  const LinearOperator<LatticeFermion>*
  UnprecWilsonFermAct::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<LatticeFermion>(linOp(state));
  }


  //! Computes the derivative of the fermionic action respect to the link field
  /*!
   *         |  dS      dS_f
   * ds_u -- | ----   + -----   ( Write )
   *         |  dU       dU
   *
   * psi -- [1./(M_dag*M)]*chi_  ( read ) 
   *
   * \param ds_u     result      ( Write )
   * \param state    gauge field ( Read )
   * \param psi      solution to linear system ( Read )
   */

  void
  UnprecWilsonFermAct::dsdu(multi1d<LatticeColorMatrix> & ds_u,
			    Handle<const ConnectState> state,
			    const LatticeFermion& psi) const
  {
    START_CODE();

    // The phi <=> X so I define the Y field as MX 
  
    // Get at the U matrices
    const multi1d<LatticeColorMatrix>& u = state->getLinks();
  
    // Get a linear operator
    Handle<const LinearOperator<LatticeFermion> > M(linOp(state));

    // Compute MY
    LatticeFermion Y;
    (*M)(Y, psi, PLUS);

    // Usually this is Kappa. In our normalisation it is 0.5 
    // I am adding in a factor of -1 to be consistent with the sign
    // convention for the preconditioned one. (We can always take this out
    // later
    Real prefactor=-Real(0.5);

    // Two temporaries
    LatticeFermion f_tmp;
    LatticeColorMatrix u_tmp;
    for(int mu = 0; mu < Nd; mu++) { 

      // f_tmp = (1 + gamma_mu) Y 
      f_tmp = Gamma(1<<mu)*Y;
      f_tmp += Y;

      //   trace_spin ( ( 1 + gamma_mu ) Y_x+mu X^{dag}_x )
//      u_tmp = traceSpin(outerProduct(shift(f_tmp, FORWARD, mu),psi));
      LatticeFermion foo = shift(f_tmp, FORWARD, mu);
      u_tmp = traceSpin(outerProduct(foo,psi));

      // f_tmp = -(1 -gamma_mu) X
      f_tmp = Gamma(1<<mu)*psi;
      f_tmp -= psi;

      //  +trace_spin( ( 1 - gamma_mu) X_x+mu Y^{dag}_x)
//      u_tmp -= traceSpin(outerProduct(shift(f_tmp, FORWARD, mu),Y));
      foo = shift(f_tmp, FORWARD, mu);
      u_tmp -= traceSpin(outerProduct(foo,Y));
    
      // accumulate with prefactor
      ds_u[mu] = prefactor*( u[mu]*u_tmp );
    }

     
    
    END_CODE();
  }

  void
  UnprecWilsonFermAct::dsdu2(multi1d<LatticeColorMatrix>& ds_u,
			    Handle<const ConnectState> state,
			    const LatticeFermion& psi) const
  {
    START_CODE();

    // The phi <=> X so I define the Y field as MX 
  
    // Get at the U matrices
    const multi1d<LatticeColorMatrix>& u = state->getLinks();
  
    // Get a linear operator
    Handle<const LinearOperator<LatticeFermion> > M(linOp(state));

    // Compute MY
    LatticeFermion Y;
    (*M)(Y, psi, PLUS);

    Real prefactor=-Real(0.5);

    // Derivative Matrix
    for(int mu=0; mu < Nd; mu++) { 

      // Create a linop for dM/dU_mu
      Handle<const LinearOperator<LatticeFermion> > dMdU(ldMdU(state,mu));

      LatticeColorMatrix u_tmp;
      LatticeFermion tmp;

      // tmp = (dM^{dagger}/dU_mu) Y
      (*dMdU)(tmp, Y, MINUS);

      // u_tmp += Tr_spin (dM^{dagger}/dU_mu) Y X^{dagger} 
      //        = traceSpin(outerProduct(tmp, X))
      u_tmp = traceSpin(outerProduct(tmp, psi));

    
      // tmp = (dM/dU_mu) X = (dM/dU_w) phi 
      (*dMdU)(tmp, psi, PLUS);

      // u_tmp = Tr_spin (dM/dU_mu) X Y^{dagger} = traceSpin(outerProduct(tmp, Y))
      u_tmp += traceSpin(outerProduct(tmp, Y));

    
      // Multiply in u_[mu]
      ds_u[mu] = u[mu]*u_tmp;

      
      

    }
    
    END_CODE();
  }



}
