// -*- C++ -*-
// $Id: hamilton.h,v 3.0 2006-04-03 04:59:10 edwards Exp $
/*! \file
 *  \brief Hamiltonian systems
 */

#ifndef __hamilton_h__
#define __hamilton_h__

#include "chromabase.h"
#include "gaugebc.h"
#include "gaugeact.h"
#include "fermact.h"


//----------------------------------------------------------------------------
//! Abstract Hamiltonian system
/*!
 * \ingroup molecdyn
 *
 * This is a crude first attempt at a Hamiltonian system.
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
 */

// OK I am a novice here, but all our fermion actions 
// are templated whereas our gauge actions are not
// so I am tempted to do something like
// Here GA is Gauge Action
// FT is the FermionType
// FA is the fermion action and it is templated on FT
template<typename GA, typename FT, class FA > 
class AbstractHamSys
{
public:
  // Derived class constructor should take params like Npf

  //! Virtual destructor
  virtual ~AbstractHamSys() {}

  //! Generate new momenta and fermionic (fixed) fields
  /*! Default version */
  virtual void newFields() 
    {
      newFerm();   // new pseudofermions
      newMom();    // new momenta
    }

  //! Copy the non-field portion of another AbstractHamSys
  virtual void copyNonFields(const AbstractHamSys& hs) = 0;

  //! Copy the field portion of another AbstractHamSys
  /*! Default version */
  virtual void copyFields(const AbstractHamSys& hs)
    {
      copyCoord(hs);
      copyMom(hs);
      copyFerm(hs);
    }

  //! These are used above, so provide virtual functions for them.
  virtual void copyCoord(const AbstractHamSys& hs) = 0;
  virtual void copyMom(const AbstractHamSys& hs) = 0;
  virtual void copyFerm(const AbstractHamSys& hs) = 0;

  //! Compute dS/dU
  virtual void dsdu(multi1d<LatticeColorMatrix>& ds_u) = 0;

  //! Refresh the momenta
  virtual void newMom() const = 0;

protected:
  //! Form the fermionic noise
  /*! 
   * Assuming here some kind of opaque state .
   * Variables like n_zero in overlap should be in the state
   * What to do with NCG?? int& n_congrd, could return it
   */
  virtual void newFerm() = 0; 
};


//! Abstract Hamiltonian system
/*!
 * \ingroup molecdyn
 *
 * This class is more specific to the case of having fermions
 * and gauge fields as coordinates.
 */

template<typename GA, typename FT,  class FA> 
class HamSys : public AbstractHamSys<GA, FT, FA>
{
public:
  // Derived class constructor should take params like Npf

  //! Virtual destructor
  virtual ~HamSys() {}

  //! Refresh the momenta
  /*! Default version */
  virtual void newMom()
    {
      /* Generate new momenta */
      for(int mu = 0; mu < Nd; ++mu)
      {
	gaussian(getMom()[mu]);

	// Destroy Linear Congruential correlators
	Real foo;
	random(foo);

	taproj(getMom()[mu]);
      }

      // Zero out the momenta on boundaries where the gauge fields
      // are held fixed.
      getGaugeAct().getGaugeBC().zero(getMom());
    }

  //! Copy the coordinate portion of another HamSys
  /*! Default version */
  virtual void copyCoord(const HamSys& hs)
    {
      getU() = hs.getU();
    }

  //! Copy the momenta portion of another HamSys
  /*! Default version */
  virtual void copyMom(const HamSys& hs)
    {
      getMom() = hs.getMom();
    }

  //! Copy the fermion portion of another HamSys
  /*! Default version */
  virtual void copyFerm(const HamSys& hs)
    {
      getPF() = hs.getPF();
      getPsi() = hs.getPsi();
    }

protected:
  //! The number of pseudo-fermions in use
  /*! This is a convenience function for getting the multi1d PF fields and using size() */
  virtual int numPF() const = 0;

  //! The gauge action in use.
  /*! The return type could instead be GA. */
  virtual const GA& getGaugeAct() const = 0;

  //! The fermion action in use
  /*! 
   * Maybe should be a precond action? Don't think so here.
   * Actually, maybe need appropriate derived FermAct to get precond. matrices.
   *
   * The return type could instead be FA .
   */
  virtual const FA& getFermAct() const = 0;

protected:
  // Somehow this stuff should be less public - maybe in a derived class
  //! Access momenta
  virtual const multi1d<LatticeColorMatrix>& getMom() const = 0;
  virtual multi1d<LatticeColorMatrix>& getMom() = 0;

  //! Access pseudofermions
  virtual const multi1d<FT>& getPF() const = 0;
  virtual multi1d<FT>& getPF() = 0;

  //! Access pseudofermions & solution vectors
  virtual const multi1d<FT>& getPsi() const = 0;
  virtual multi1d<FT>& getPsi() = 0;

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
template<typename GA, typename FT,  class FA> 
class ExactHamSys : public HamSys<GA, FT, FA >
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
  //! Form the fermionic noise
  virtual void newFerm()
    {
      // The connection state (of the gauge field) includes any fermion BC
      Handle<const ConnectState> state(getFermAct().createState(getU()));

      // Get the linop
      /*! WARNING: need a generic template param here */
      Handle<const LinearOperator<FT> > A(getFermAct().linOp(state));

      // Generate the pseudofermions
      for(int i=0; i < numPF(); ++i)
      {
	// Initialize psi
	FT psi = zero;
  
	// Fill eta with random gaussian noise: < eta_dagger * eta > = 1
	FT eta = zero;
	gaussian(eta, A.subset());
  
	// Zero out any potential constant boundaries
	// BJ: DO WE NEED THIS? FermAct doesnt have a zero() method
	// RGE: fixed
	getFermAct().getFermBC().zero(eta);  // WARNING: WHAT ABOUT SUBSET HERE???

	// chi = M_dag*eta
	FT chi;
	A(chi, eta, MINUS);

	// Insert into the HamSys
	getPF()[i] = chi;
      }
    }

  //! Compute new Psi fields
  /*! 
   * This may be called lazily and the result cached, 
   * e.g., only if the U fields change 
   */
  virtual void newPsi()
    {
      // The connection state (of the gauge field) includes any fermion BC
      Handle<const ConnectState> state(getFermAct().createState(getU()));

      // Get the linop
      /*! WARNING: need a generic template param here */
      Handle<const LinearOperator<FT> > A(getFermAct().linOp(state));

      // Solve the linear system for the psi fields
      for(int i=0; i < numPF(); ++i)
      {
	// Initialize psi
	FT psi = zero;

	// psi = (1/(M_dag*M))*chi
	InvCG2(*A, getPF()[i], psi, getRsdCG(), getMaxCG(), n_count);

	// Out of paranoia, zero out the boundaries here
	/* NOTE: This may not be needed */
	// BJ: No zero method in FermAct
	//
	//getFermAct().zero(psi);  // WARNING: WHAT ABOUT SUBSET HERE???

	// Insert into the HamSys
	getPsi()[i] = psi;
      }
    }

  //! Get the residual
  virtual Real getRsdCG(void) const = 0;

  //! Get the max number of iterations
  virtual int  getMaxCG(void) const = 0;

  //! Measure the potential energy
  virtual Double mesPE(void) = 0;

  //! Measure the kinetic energy
  virtual Double mesKE(void) = 0;

  //! Measure the fermionic energy
  virtual Double mesFE(void) = 0;
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

template<typename GA, typename FT, class FA> 
class SingleExactHamSys : public ExactHamSys<GA, FT, FA>
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


protected:
  //! Compute dS_g/dU
  /*! This routine should include all scale factors */
  virtual void dsdu_g(multi1d<LatticeColorMatrix>& ds_u)
    {
      // The connection state (of the gauge field) includes any gauge BC
      Handle<const ConnectState> state(getGaugeAct().createState(getU()));

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

      ds_u = 0;
      for(int i=0; i < numPF(); ++i)
      {
	getFermAct().dsdu(ds_u, state, getPsi()[i]);

	for(int mu=0; mu < Nd; ++mu)
	  ds_u[mu] += tmp[mu];
      }

//    for(int mu=0; mu < Nd; ++mu)
//      ds_u[mu] *= factors;     // possibly put in some ferm factors??
    }

};

#if 0
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

  //! Virtual destructor
  virtual ~RatExactHamSys() {}

  // More blah blah blah
  // Now the dsdu_f is on the rat. approx. over the fermion action, which
  // at some point needs the fermion dsdu_f.
};
#endif


//---------------------------------------------------------------------------

//! Type of trajectory termination
enum AlgETrj_t {FIXED_LENGTH, EXPONENTIAL_LENGTH};


//! Abstract molecular dynamics single trajectory integrator
/*!
 * \ingroup molecdyn
 *
 * This is a crude first attempt at a molecular dynamics integrator
 * Lots more thought needed here.
 */

template<typename GA, typename FT, typename FA, template<typename,typename,typename> class HS>
class HybInt
{
public:
  //! Virtual destructor
  virtual ~HybInt() {}

  //! Do an integration
  virtual void operator()(HS<GA,FT,FA>& ham) const = 0;

  //! Leaps P forward step eps1 & eps2
  /*! This default version is appropriate for all HamSys types */
  virtual void leapPMX(HS<GA,FT,FA>& ham, 
		       const Real& eps1, const Real& eps2,
		       bool EndP) const 
  {
      if (EndP)    // the szin version has many things here
	return;

      Real eps12 = eps1 + eps2;
      LeapP(ham, eps12);
  }

  //! Leap a step in the momenta
  /*! Default version */
  virtual void leapP(HS<GA,FT,FA>& ham, const Real& eps) const
  {
      multi1d<LatticeColorMatrix> ds_u(Nd);
      /* Derivative of action (includes gauge and ferm !!)
       * This must include knowledge of say poly approx.
       * Also, ham should be ``Exact'' or ``Inexact'' for Phi or R alg's
       */
      int niter = 0;
      ham.dsdu(ds_u);

      /* p_mom = p_mom - eps*[Ds/Du + DS_f/Du] */
      for(int mu=0; mu<Nd; ++mu)
	ham.getMom()[mu] -= eps * ds_u[mu];
      
      for(int mu=0; mu<Nd; ++mu)
	taproj(ham.getMom()[mu]);

//      return niter;
  }

  //! Leap a step in the gauge fields
  /*! Default version */
  virtual void leapU(HS<GA,FA>& ham, const Real& eps) const
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

      // Apply the BC to the gauge fields
      /* NOTE: this may not be necessary if the reunit  */
      getGaugeAct().modify(u);  // Something like this is used in SZIN
  }

  //! Return the step size
  virtual Real getStepSize() const = 0;

  //! Return the expected trajectory length 
  /*! In same units as step-size, e.g. NOT multiples of step-size */
  virtual Real getTrjLen() const = 0;

  //! Type of termination for the end of trajectory
  virtual AlgETrj_t AlgETrjType() const = 0;    // Specify trajectory length algorithm

  //! Test for end of a trajectory
  bool endOfTrj(const Real& t) const
  {
      Real r;
      Real dt = getStepSize();
  
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
	QDPIO::cerr << "endOfTrj: unknown algorithm termination" << endl;
	QDP_abort(1);
      }
    }
};



#if 0 
//! Abstract leap-frog molecular dynamics single trajectory integrator
/*!
 * \ingroup molecdyn
 *
 *  NOTE: this class is not derived from exact integrator. 
 */
template<class GA, class FA, template<class,int> class HS>
class LeapFrogHybInt : public HybInt<GA,FA,HS>
{
public:
  //! Virtual destructor
  virtual ~LeapFrogHybInt() {}

  //! Do an integration
  virtual void operator()(HS<GA,FA>& ham) const
    {
      bool EndP = false;
      Real dt = getStepSize();
      Real dtH = 0.5*dt;
      Real t = 0;

      int niter = 0;
      LeapP(ham, dtH);   // first half step

      while (! EndP)
      {
	t += dt;

	LeapU(ham, dt);

	// NOTE: SZIN uses an interpolation of the Psi field here

	EndP = endOfTrj(t);
	LeapPMX(han, dtH, dtH, EndP);  // full step

	// NOTE: no monitoring here
      }

//      return niter;
    }
};





//! Abstract Campostrini molecular dynamics single trajectory integrator
/*!
 * \ingroup molecdyn
 *
 *  NOTE: this class is not derived from exact integrator. 
 */
template<class GA, class FA, template<class,int> class HS>
class CampHybInt : public HybInt<GA,FA,HS>
{
public:
  //! Virtual destructor
  virtual ~CampHybInt() {}

  //! Do an integration
  virtual void operator()(HS<GA,FA>& ham) const
    {
      bool EndP = false;
      Real dt = getStepSize();
      Real t = 0;
      int niter = 0;

      Real sigma    = pow(2.0, 1.0/3.0);
      Real eps      = dt / (2 - sigma);
      Real epsH     = 0.5 * eps;
      Real epsuBac  = -sigma * eps;
      Real epsSigH  = -sigma * epsH;
      Real sigmaInv = 1.0 / sigma;

      LeapP(ham, epsH);   // first half step

      while (! EndP)
      {
	t += dt;
	int wgl_cnt = 0;

	// First step of wiggle
	LeapU(ham, eps);

	// WARNING: What about interpolation here???
	// Interpol(psi, old_psi, Real(-1), Npf);

	LeapPMX(ham, epsH, epsSigH, false);

	// NOTE: no monitoring here

	// Second step of wiggle
	LeapU(ham, epsUBac);

	// WARNING: What about interpolation here???
	// interpol(psi, old_psi, sigma, Npf);

	LeapPMX(ham, epsSigH, epsH, false);

	// NOTE: no monitoring here

	// Third step of wiggle
	LeapU(ham, eps);

	// WARNING: What about interpolation here???
	// interpol(psi, old_psi, sigmainv, Npf);

	EndP = endOfTrj(t);
	LeapPMX(ham, epsH, epsH, EndP);

	// NOTE: no monitoring here
      }

//      return niter;
    }
};
#endif 


//! Abstract molecular dynamics single trajectory integrator for exact actions
/*!
 * \ingroup molecdyn
 *
 * The inheritance here is a nuisance, I want leapPMX generic to leapfrog
 * and Campostrini
 */
template<class GA, class FA, template<class,int> class HS>
class ExactLeapFrogHybInt : public LeapFrogHybInt<GA,FA,HS>
{
public:
  //! Virtual destructor
  virtual ~ExactHybInt() {}

  //! Monitor the integration?
  /*! Default version */
  virtual bool monitordH() const {return false;}

  //! Leaps P forward step eps1 & eps2 with optional monitoring of dH
  /*! This default version is appropriate for simple/single HamSys types */
  virtual void leapPMX(HS<GA,FA>& ham, 
		       const Real& eps1, const Real& eps2,
		       bool EndP) const
    {
//      if (EndP)    // the szin version has many things here
//	return;

      int niter = 0;
      Real eps12 = eps1 + eps2;

      /*
       * The getPsi could/should be lazy, but I will help it along here.
       * Maybe this should go away since the idea of a psi is an optimization.
       * We only need inverses when they are requested, namely for energies
       * and force calcs.
       * However, SZIN's leappmx exposes this to use for a trick, namely
       * to watch the acceptance on a single step change as the resid. is
       * decreased. This is too much an optimization to be used for a 
       * first default version.
       */
      if ( ! monitordH() && ! EndP )
      {
	LeapP(ham, eps12);
      }
      else
      {
	LeapP(ham, eps1);
	if ( monitordH() )
	{
//	  ham.mesE(ke,pe,fe);   // problem - wants to return E's in middle of traj.
	}
	LeapP(ham, eps2);
      }
      
//      return niter;
    }
};


//----------------------------------------------------------------------------------
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
  virtual bool operator()(HS<GA,FA>& hmc, HS<GA,FA>& hmd) const
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

#endif
