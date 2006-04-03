// -*- C++ -*-
// $Id: field_state.h,v 3.0 2006-04-03 04:59:07 edwards Exp $

/*! \file
 * \brief Field state
 *
 * Abstract classes to define the states of fields
 * This is basically a collection class...
 *
 * Templated to allow generic situations
 */

#ifndef field_state_h
#define field_state_h

#include "chromabase.h"

namespace Chroma 
{
  //! Abstract field state
  /*! @ingroup molecdyn
   *
   * An abstract field state. The templates P and Q are the 
   * types of the  coordinates and momenta 
   */
  template <typename P, typename Q>
  class AbsFieldState 
  {
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
  
  /*
    template <typename P, typename Q, typename Phi>
    class AbsPFFieldState : public AbsFieldState<P,Q> {
    public:

    //! Virtual destructor
    virtual ~AbsPFFieldState<P,Q,Phi>() {}

    //! Clone the state -- this will be needed for the virtual 
    // copy ideas for derived classes
    virtual AbsPFFieldState<P,Q,Phi>* clone(void) const = 0;

    //! Accessors
    virtual const P& getP(void) const = 0;
    virtual const Q& getQ(void) const = 0;

    virtual const int numPhi(void) const =0;
    virtual const Phi& getPhi(int i) const = 0;

    //! Mutators
    virtual P& getP(void) = 0;
    virtual Q& getQ(void) = 0;
    virtual Phi& getPhi(int i) = 0;

    // And that is all we can do here
    };
  */

  //! Pure gauge field state
  /*! @ingroup molecdyn
   *
   * A pure Gauge field state type field state. 
   */
  class GaugeFieldState : public AbsFieldState<multi1d<LatticeColorMatrix>,
			                       multi1d<LatticeColorMatrix> >
  {
  public: 
      
    // Constructor
    GaugeFieldState(const multi1d<LatticeColorMatrix>& p_,
		    const multi1d<LatticeColorMatrix>& q_) {
      p.resize(Nd);
      q.resize(Nd);
      for(int mu=0; mu < Nd; mu++) { 
	p[mu] = p_[mu];
	q[mu] = q_[mu];
      }
    }
      
    // Copy Constructor
    GaugeFieldState(const GaugeFieldState& s)  {
      p.resize(Nd);
      q.resize(Nd);
      for(int mu=0; mu < Nd; mu++) { 
	p[mu] = s.p[mu];
	q[mu] = s.q[mu];
      }
	
    }
      
    // Destructor
    ~GaugeFieldState() {};

    // Clone function -- covariant return type
    GaugeFieldState* clone(void) const { 
      return new GaugeFieldState(*this);
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
  
};  // End namespace Chroma

#endif
