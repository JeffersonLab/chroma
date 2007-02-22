// -*- C++ -*-
// $Id: pg_gaugeact.h,v 3.2 2007-02-22 21:11:48 bjoo Exp $
/*! \file
 *  \brief Parallelgram gauge action
 */

#ifndef __pg_gaugeact_h__
#define __pg_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"

namespace Chroma
{

  /*! @ingroup gaugeacts */
  namespace PgGaugeActEnv 
  { 
    extern const string name;
    bool registerAll();
  }

  // Parameter structure
  /*! @ingroup gaugeacts */
  struct PgGaugeActParams {
    // Base Constructor
    PgGaugeActParams();
    
    // Read params from some root path
    PgGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real coeff;  
  };
  
  /*! @ingroup gaugeacts */
  void read(XMLReader& xml, const string& path, PgGaugeActParams& param);
  

  //! Parallelogram gauge action
  /*! \ingroup gaugeacts
   *
   * The standard parallelogram gauge action
   */

  class PgGaugeAct : public LinearGaugeAction
  {
  public:
    //! General CreateGaugeState
    PgGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
	       const Real& coeff_) : 
      cgs(cgs_), coeff(coeff_) {}

    //! Read coeff from a param struct
    PgGaugeAct(Handle< CreateGaugeState<P,Q> > cgs_, 
	       const PgGaugeActParams& p) :
      cgs(cgs_), coeff(p.coeff) {}

    //! Is anisotropy used?
    bool anisoP() const {return false;}

    //! Anisotropy factor
    const Real anisoFactor() const {return Real(1);}

    //! Anisotropic direction
    int tDir() const {return Nd-1;}

    //! Return the set on which the gauge action is defined
    /*! Defined on the even-off (red/black) set */
    const Set& getSet() const {return rb;}

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void staple(LatticeColorMatrix& result,
		const Handle< GaugeState<P,Q> >& state,
		int mu, int cb) const;

    //! Compute dS/dU
    void deriv(multi1d<LatticeColorMatrix>& result,
	       const Handle< GaugeState<P,Q> >& state) const;

    //! Compute the actions
    Double S(const Handle< GaugeState<P,Q> >& state) const;

    //! Produce a gauge create state object
    const CreateGaugeState<P,Q>& getCreateState() const {return *cgs;}

    //! Destructor is automatic
    ~PgGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getCoeff(void) const { return coeff; }

  protected:
    //! Partial constructor
    PgGaugeAct() {}
    //! Hide assignment
    void operator=(const PgGaugeAct& a) {}

  private:
    Handle< CreateGaugeState<P,Q> >  cgs;  // create gauge state
    Real coeff;              // The coupling coefficient

  };

};


#endif
