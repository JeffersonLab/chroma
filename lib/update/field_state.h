#ifndef field_state_h
#define field_state_h

/*! Abstract classes to define the states of fields
  This is basically a collection class...

  Templated to allow generic situations
*/
#include "chromabase.h"
using namespace QDP;

/*! An abstract field state. The templates P and Q are the 
  types of the  coordinates and momenta */
template <typename P, typename Q>
class AbsFieldState {
 public:

  //! Virtual destructor
  virtual ~AbsFieldState<P,Q>() {}

  //! Clone the state -- this will be needed for the virtual 
  // copy ideas for derived classes
  virtual AbsFieldState<P,Q>* clone(void) const = 0;

  //! Accessors
  virtual const P& getP(void) const = 0;
  virtual const Q& getQ(void) const = 0;

  //! Mutators
  virtual P& getP(void) = 0;
  virtual Q& getQ(void) = 0;

  // And that is all we can do here
};

/*! An abstract fermionic field state. This does the same as 
  AbsFieldState, but now allows for pseudofermions. Note I am not
  allowing refresh as that would need to know the action.

  Since it is awkward to fill the copy() method of a non fermionic
  AbsFieldState I don't derive from it. This is kinda like the Lattice
  Fermion and LatticeFermion<multi1d<>> construction 
*/
template <typename P, typename Q, typename Phi>
class AbsPFFieldState {
 public:

  //! Virtual destructor
  virtual ~AbsPFFieldState<P,Q,Phi>() {}

  // Copy for later virtual copying 
  virtual AbsPFFieldState<P,Q,Phi>* clone(void) const = 0;
 
  // Accessors
  virtual const P&             getP(void)  const = 0;
  virtual const Q&             getQ(void)  const = 0;
  virtual const multi1d<Phi>&  getPF(void) const = 0;

  // Mutators
  virtual P&                   getP(void) = 0;
  virtual Q&                   getQ(void) = 0;
  virtual multi1d<Phi>&        getPF(void) = 0;
};

/*! A pure Gauge field state
  type field state. */
class PureGaugeFieldState : public AbsFieldState<multi1d<LatticeColorMatrix>,                                                    multi1d<LatticeColorMatrix> >
{
 public: 

  // Constructor
  PureGaugeFieldState(multi1d<LatticeColorMatrix>& p_,
		      multi1d<LatticeColorMatrix>& q_) : p(p_), q(q_) {}

  // Copy Constructor
  PureGaugeFieldState(const PureGaugeFieldState& s) : p(s.p), q(s.q) {}

  // Destructor
  ~PureGaugeFieldState() {};

  // Clone function -- covariant return type
  PureGaugeFieldState* clone(void) const { 
    return new PureGaugeFieldState(*this);
  }

  // Accessors
  const multi1d<LatticeColorMatrix>& getP(void) const { return p; }
  const multi1d<LatticeColorMatrix>& getQ(void) const { return q; }

  // Mutators
  multi1d<LatticeColorMatrix>& getP(void)  { return p; }
  multi1d<LatticeColorMatrix>& getQ(void)  { return q; }
  
 private:
  multi1d<LatticeColorMatrix> p;
  multi1d<LatticeColorMatrix> q;
};


/*! Try a gauge,latticeFermion type field state. */
class GaugeFermFieldState : public AbsPFFieldState<multi1d<LatticeColorMatrix>,                                                    multi1d<LatticeColorMatrix>,
                                                   LatticeFermion>
{
 public: 

  // Constructor
  GaugeFermFieldState(multi1d<LatticeColorMatrix>& p_,
		      multi1d<LatticeColorMatrix>& q_,
		      multi1d<LatticeFermion>& pf_) : p(p_), q(q_), pf(pf_) {}

  // Copy Constructor
  GaugeFermFieldState(const GaugeFermFieldState& s) : p(s.p), q(s.q), pf(s.pf) {}

  // Destructor
  ~GaugeFermFieldState() {};

  // Clone function -- covariant return type
  GaugeFermFieldState* clone(void) const { 
    return new GaugeFermFieldState(*this);
  }

  // Accessors
  const multi1d<LatticeColorMatrix>& getP(void) const { return p; }
  const multi1d<LatticeColorMatrix>& getQ(void) const { return q; }
  const multi1d<LatticeFermion>&     getPF(void) const { return pf; }

  // Mutators
  multi1d<LatticeColorMatrix>& getP(void)  { return p; }
  multi1d<LatticeColorMatrix>& getQ(void)  { return q; }
  multi1d<LatticeFermion>& getPF(void) { return pf; }
  
 private:
  multi1d<LatticeColorMatrix> p;
  multi1d<LatticeColorMatrix> q;
  multi1d<LatticeFermion> pf;
};
#endif
