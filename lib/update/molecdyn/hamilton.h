// -*- C++ -*-
// $Id: hamilton.h,v 1.5 2003-12-31 23:59:28 edwards Exp $
/*! \file
 *  \brief Hamiltonian systems
 */

#ifndef __hamilton_h__
#define __hamilton_h__

#include "fermact.h"

using namespace QDP;

//! Abstract Hamiltonian system
/*!
 * \ingroup molecdyn
 *
 * This is a really crude first attempt at a Hamiltonian system.
 * Lots more thought needed here.
 *
 * Maybe these functions should be global? No, I think there is a case
 * to be made for having an "MD Integrator" class with funcs.
 * An HMDIntegrator should take a HamSys template param'd by action types.
 * I suspect I should have a typelist of actions, (e.g., can have
 * more than 1 ferm action). The HamSys constructor would then take
 * instances of those actions. 
 *
 * There should be some general notion of "Coordinates" and "Conjugate Momenta".
 * How pseudo-ferm cleanly fit in that picture in this code is not clear
 * (to me, at least).
 *
 * NOTE: dsdu is found in the fermion (or gauge) action
 */

template<class GA, class FA> 
class HamSys
{
public:
  // Derived class constructor should take params like Npf

  //! Virtual destructor
  virtual ~HamSys() {}

  //! Generate new momenta and fermionic (fixed) fields
  /*! Default version */
  virtual void newFields()
    {
      newFerm();   // new pseudofermions
      newMom();    // new momenta
    }

  //! Copy the non-field portion of another HamSys
  virtual void copyNonFields(const HamSys& hs) = 0;

  //! Copy the field portion of another HamSys
  /*! Default version */
  virtual void copyFields(const HamSys& hs)
    {
      copyCoord(hs);
      copyMom(hs);
      copyFerm(hs);
    }

  //! Compute dS/dU
  virtual void dsdu(multi1d<LatticeColorMatrix>& ds_u) = 0;

protected:
  //! Refresh the momenta
  virtual void newMom() = 0;

  //! Form the fermionic noise
  /*! 
   * Assuming here some kind of opaque state .
   * Variables like n_zero in overlap should be in the state
   * What to do with NCG?? int& n_congrd, could return it
   */
  virtual void newFerm() = 0; 

#if 0
// Probably shouldn't be here
  //! Compute dS/dU
  virtual void dsdu_g(multi1d<LatticeColorMatrix>& ds_u,
		    multi1d<LatticeColorMatrix>& u) const = 0;

  //! Compute dS_f/dU
  virtual void dsdu_f(multi1d<LatticeColorMatrix>& result,
		      multi1d<LatticeColorMatrix>& u,
		      const multi1d<LatticeFermion>& psi) const = 0;
#endif

public:
  // Somehow this stuff should be less public - maybe in a derived class
  //! Access momenta
  virtual const multi1d<LatticeColorMatrix>& getMom() const = 0;
  virtual multi1d<LatticeColorMatrix>& getMom() = 0;

  //! Access pseudofermions
  virtual const multi1d<LatticeFermion>& getPF() const = 0;
  virtual multi1d<LatticeFermion>& getPF() = 0;

  //! Access pseudofermions & solution vectors
  virtual const multi1d<LatticeFermion>& getPsi() const = 0;
  virtual multi1d<LatticeFermion>& getPsi() = 0;

  //! Access u fields
  virtual const multi1d<LatticeColorMatrix>& getU() const = 0;
  virtual multi1d<LatticeColorMatrix>& getU() = 0;
};


//! Abstract HamSys with that supports exact calc. of energy
/*!
 * \ingroup molecdyn
 *
 *  Suitable for use in a "Phi" HMD/HMC algorithm.
 *
 */

template<class GA, class FA> 
class ExactHamSys : public HamSys<GA,FA>
{
public:
  // Derived class constructor should take params like Npf

  //! Virtual destructor
  virtual ~ExactHamSys() {}

  //! Measure the total energy of the system
  /*! 
   * The energy measurements use the internally held coord/mom. The idea
   * used here is that it only makes sense to ask for these energies when all
   * the fields are at the same step.
   *
   * The energies are held in some internal state. They may not have to
   * be recomputed if fields have not changed since a last calculation.
   */
  virtual void mesE(Double& ke, Double& pe, Double& fe)
    {ke = mesKE(); pe = mesPE(); fe = mesFE();}

protected:
  //! Measure the potential energy
  virtual Double mesPE() = 0;

  //! Measure the kinetic energy
  virtual Double mesKE() = 0;

  //! Measure the fermionic energy
  virtual Double mesFE() = 0;
};


//! Abstract HamSys with simple ferm. form and exact calc.
/*!
 * \ingroup molecdyn
 *
 *  Ham. Sys. of form   H = S_G + \chi^dagger [M(U)]^{-1} \chi
 *  Energy can be computed exactly.
 *  The "single" (okay, come up with a better name...) is this simple
 *  ferm. form, as opposed to poly. approx. or rat. approx.
 */

template<class GA, class FA> 
class SingleExactHamSys : public ExactHamSys<GA,FA>
{
public:
  // Derived class constructor should take params like Npf

  //! Virtual destructor
  virtual ~SingleExactHamSys() {}

  //! Compute dS/dU
  virtual void dsdu(multi1d<LatticeColorMatrix>& ds_u)
    {
      dsdu_g(ds_u);   // E.g., this must include BetaMD  !!

      multi1d<LatticeColorMatrix> ds_u_f(Nd);
      dsdu_f(ds_u_f); // Here is where different ferm. forms. exploited

      // Add on the appropriately scaled gauge and fermion dsdu
      for(int mu=0; mu < Nd; ++mu)
	ds_u[mu] += ds_u_f[mu];

      gauge().zero(ds_u);   // Zero out any potential constant boundaries
    }

  //! Compute dS_g/dU
  /*! This routine should include all scale factors */
  virtual void dsdu_g(multi1d<LatticeColorMatrix>& ds_u)
    {
      // The connection state (of the gauge field) includes any gauge BC
      Handle<const ConnectState> state(getFermAct().createState(getU()));

      getGaugeAct().dsdu(ds_u, state);   // E.g., this must include BetaMD  !!
    }

  //! Compute dS_f/dU
  /*! 
   * Here is where different ferm. forms. exploited. 
   * This routine should include all scale factors.
   */
  virtual void dsdu_f(multi1d<LatticeColorMatrix>& ds_u)
    {
      // The connection state (of the gauge field) includes any fermion BC
      Handle<const ConnectState> state(getFermAct().createState(getU()));

      // Evaluate dsdu_f with this fermact on this state and psi
      multi1d<LatticeColorMatrix> tmp(Nd);

      getFermAct().dsdu(ds_u, state, getPsi()[0]);
      for(int i=1; i < numPF(); ++i)    // unroll the loop
	for(int mu=0; mu < Nd; ++mu)
	  ds_u[mu] += tmp[mu];

//    for(int mu=0; mu < Nd; ++mu)
//      ds_u[mu] *= factors;     // possibly put in some ferm factors??
    }

  //! The gauge action in use.
  /*! The return type could instead be GA. */
  virtual const GaugeAction& getGaugeAct() const = 0;

  //! The fermion action in use
  /*! 
   * Maybe should be a precond action? Don't think so here.
   * Actually, maybe need appropriate derived FermAct to get precond. matrices.
   *
   * The return type could instead be FA .
   */
  virtual const FermionAction<LatticeFermion>& getFermAct() const = 0;
};


//! Abstract HamSys with rational approx to ferm form and exact calc.
/*!
 * \ingroup molecdyn
 *
 *  Ham. Sys. of form   H = S_G + \chi^dagger [P(U)/Q(U)]^{-1} \chi
 *  Energy can be computed exactly.
 *  The "single" (okay, come up with a better name...) is this simple
 *  ferm. form, as opposed to poly. approx. or rat. approx.
 *
 *  RatExact is a cool name, but maybe ExactRat is better?
 */

template<class GA, class FA> 
class RatExactHamSys : public ExactHamSys<GA,FA>
{
public:
  // Derived class constructor should take params like Npf

  // More blah blah blah
  // Now the dsdu_f is on the rat. approx. over the fermion action, which
  // at some point needs the fermion dsdu_f.

  //! Virtual destructor
  virtual ~RatExactHamSys() {}
};



//! Abstract molecular dynamics single trajectory integrator
/*!
 * \ingroup molecdyn
 *
 * This is a really crude first attempt at a molecular dynamics integrator
 * Lots more thought needed here.
 *
 * An HMDIntegrator should take a HamSys template param'd by action types.
 * I suspect I should have a typelist of actions, (e.g., can have
 * more than 1 ferm action). The HamSys constructor would then take
 * instances of those actions. 
 *
 * There should be some general notion of "Coordinates" and "Conjugate Momenta".
 * How pseudo-ferm cleanly fit in that picture in this code is not clear
 * (to me, at least).
 */

template<class GA, class FA, template<class,int> class HS>
class HybInt
{
public:
  //! Do an integration
  virtual void operator()(HS<GA,FA>& ham) const = 0;

  //! Leaps P forward step eps with optional monitor of DelH
  virtual void leapPMX(HS<GA,FA>& ham) const = 0;

  //! Leap a step in the momenta
  /*! Default version */
  virtual void leapP(HS<GA,FA>& ham, const Real& eps) const
    {
      multi1d<LatticeColorMatrix> ds_u(Nd);
      /* Derivative of action (includes gauge and ferm !!)
       * This must include knowledge of say poly approx.
       * Also, ham should be ``Exact'' or ``Inexact'' for Phi or R alg's
       */
      ham.dsdu(ds_u);

      /* p_mom = p_mom - eps*[Ds/Du + DS_f/Du] */
      for(int mu=0; mu<Nd; ++mu)
	ham.getMom()[mu] -= eps * ds_u[mu];
      
      for(int mu=0; mu<Nd; ++mu)
	taproj(ham.getMom()[mu]);
    }

  //! Leap a step in the gauge fields
  /*! Default version */
  virtual void leapU(HS<GA,FA>& ham) const
    {
      LatticeColorMatrix tmp_1;
      LatticeColorMatrix tmp_2;

      for(int mu=0; mu < Nd; ++mu)
      {
	/* tmp_1 = eps*p_mom */
	tmp_1 = p_mom[mu] * eps;

	/* tmp_1 = e^(tmp_1) */
	/* Method of exp buried here */
	expmat(tmp_1);

	/* tmp_2 = u_old * tmp_1 = u_old * exp[ eps * p_mom ] */
	tmp_2 = tmp_1 * u[mu];
  
	/* Fill new u and reunitarize ---- */
	u[mu] = tmp_2;
	int numbad;
	reunit(u[mu], REUNITARIZE_ERROR, numbad); 
      }
    }

  //! Return expected trajectory length
  virtual Real getTrjLen() const = 0;

  //! Type of termination for the trajectory
  virtual AlgETrj_t AlgETrjType() const = 0;    // Specify trajectory length algorithm

  //! Test for end of a trajectory
  bool EndOfTrj(const Real& t) const
    {
      Real r;
  
      switch (AlgETrjType())
      {
      case FIXED_LENGTH:
	EndP = ( toBool(t > 0.5*(getTrjLen()+dt)) ) ? true : false;
	break;

      case EXPONENTIAL_LENGTH:
	random(r);
	EndP = ( toBool(r < dt/max(getTrjLen()-dt,dt)) ) ? true : false;
	break;

      default:
	QDP_error_exit("unknown algorithm termination", AlgETrj);
      }
    }


  //! Virtual destructor
  virtual ~HybInt() {}
};


//! Abstract single trajectory Hybrid Monte Carlo(/HMD) integrator
/*!
 * \ingroup molecdyn
 *
 *  The main actual HMC integrator, including MD integrator
 */

template<class GA, class FA, template<class,int> class HS>
class HMCTrj
{
public:
  // Constructor should take a MD integrator as an arg

  //! Do an integration
  virtual void operator()(HS<GA,FA>& hmc, HS<GA,FA>& hmd) const
    {
      hmc.newFields();       // new mom, ferm fields
      hmd.copyFields(hmc);   // coord, mom, ferm fields
      getHybInt()(hmd); // do one HMD integration

      HS<GA,FA> tmp;               // problematic, HS may not be concrete...
      tmp.copyNonFields(hmc);
      tmp.copyFields(hmd);   // coord, mom, ferm fields
      bool acc = accept(hmc,tmp);  // accept/reject step with tmp HS
      return acc;
    }

  //! Accept/reject step
  virtual bool accept(HS<GA,FA>& hmc, HS<GA,FA>& tmp) const
    {
      Double ke_old, pe_old, fe_old;
      hmc.mesE(ke_old, pe_old, fe_old);    // old energy

      Double ke_new, pe_new, fe_new;
      tmp.mesE(ke_new, pe_old, fe_old);    // new energy

      DelKe = ke_new - ke_old; DelPe = pe_new - pe_old; DelFe = fe_new - fe_old;
      Double DelH = DelKe - DelPe + DelFe;

      // Accept/reject part
      bool acc;
      if (DelH < 0)
	acc = true;
      else
      {
	Double r;
	random(r);
	acc = (toBool(DelH < r)) ? true : false;
      }
      if (acc)
	hmc.fields() = tmp.fields();  // maybe just the coord ??

      return acc;
    }

  //! Access the MD integrator
  const HybInt<GA,FA,HS>& getHybInt() const = 0;

  //! Virtual destructor
  virtual ~HMCTrj() {}
};



#endif
