// -*- C++ -*-
// $Id: pg_gaugeact.h,v 1.1 2005-01-13 02:02:38 edwards Exp $
/*! \file
 *  \brief Parallelgram gauge action
 */

#ifndef __pg_gaugeact_h__
#define __pg_gaugeact_h__

#include "gaugeact.h"
#include "gaugebc.h"

namespace Chroma
{

  namespace PgGaugeActEnv { 
    extern const string name;
    extern const bool registered;
  }

  // Parameter structure
  struct PgGaugeActParams {
    // Base Constructor
    PgGaugeActParams();
    
    // Read params from some root path
    PgGaugeActParams(XMLReader& xml_in, const std::string& path);

    Real coeff;  
  };
  
  void read(XMLReader& xml, const string& path, PgGaugeActParams& param);
  

  //! Parallelogram gauge action
  /*! \ingroup gaugeact
   *
   * The standard parallelogram gauge action
   */

  class PgGaugeAct : public GaugeAction
  {
  public:
    //! General GaugeBC
    PgGaugeAct(Handle< GaugeBC > gbc_, 
	       const Real& coeff_) : 
      gbc(gbc_), coeff(coeff_) {}

    //! Read coeff from a param struct
    PgGaugeAct(Handle< GaugeBC > gbc_, 
	       const PgGaugeActParams& p) :
      gbc(gbc_), coeff(p.coeff) {}

    //! Copy constructor
    PgGaugeAct(const PgGaugeAct& a) : 
      gbc(a.gbc), coeff(a.coeff) {}


    //! Assignment
    PgGaugeAct& operator=(const PgGaugeAct& a)
    {gbc=a.gbc; coeff=a.coeff; return *this;}

    //! Is anisotropy used?
    bool anisoP() const {return false;}

    //! Anisotropy factor
    const Real anisoFactor() const {return Real(1);}

    //! Anisotropic direction
    int tDir() const {return Nd-1;}

    //! Return the set on which the gauge action is defined
    /*! Defined on the even-off (red/black) set */
    const OrderedSet& getSet() const {return rb;}

    //! Produce a gauge boundary condition object
    const GaugeBC& getGaugeBC() const {return *gbc;}

    //! Compute staple
    /*! Default version. Derived class should override this if needed. */
    void staple(LatticeColorMatrix& result,
		Handle<const ConnectState> state,
		int mu, int cb) const;

    //! Compute dS/dU
    void dsdu(multi1d<LatticeColorMatrix>& result,
	      const Handle<const ConnectState> state) const;

    //! Compute the actions
    Double S(const Handle<const ConnectState> state) const;

    //! Destructor is automatic
    ~PgGaugeAct() {}

    // Accessors -- non mutable members.
    const Real getCoeff(void) const { return coeff; }

  private:
    Handle< GaugeBC >  gbc;  // Gauge Boundary Condition
    Real coeff;              // The coupling coefficient

  };

};

using namespace Chroma;

#endif
