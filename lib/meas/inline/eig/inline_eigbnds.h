// -*- C++ -*-
// $Id: inline_eigbnds.h,v 1.3 2005-04-06 04:34:53 edwards Exp $

/*! \file
 * \brief Inline measurements for eigenvalue bounds
 *
 * Measure bounds of M^dag*M
 */

#ifndef __inline_eigbndsmdagm_h__
#define __inline_eigbndsmdagm_h__

#include "chromabase.h"
#include "fermact.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  namespace InlineEigBndsMdagMEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Parameter structure
  struct InlineEigBndsMdagMParams 
  {
    InlineEigBndsMdagMParams();
    InlineEigBndsMdagMParams(XMLReader& xml_in, const std::string& path);

    unsigned long frequency;
    Handle< const FermionAction<LatticeFermion> > fermact;

    //! Struct for parameters needed for a Ritz type solve
    struct RitzParams_t
    {
      Real RsdR;
      Real RsdA;
      Real RsdZero;
      bool ProjApsiP;
      int  Nmin;
      int  MaxCG;
      int  Nrenorm;
    } ritz;
  };


  //! Inline measurement of eigenvalue bounds of M^dag*M
  class InlineEigBndsMdagM : public AbsInlineMeasurement 
  {
  public:
    ~InlineEigBndsMdagM() {}
    InlineEigBndsMdagM(const InlineEigBndsMdagMParams& p) : params(p) {}
    InlineEigBndsMdagM(const InlineEigBndsMdagM& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! 4D
    void do4d(Handle< const LinearOperator<LatticeFermion> > MM,
	      unsigned long update_no,
	      XMLWriter& xml_out);

    //! 5D
    void do5d(Handle< const LinearOperator< multi1d<LatticeFermion> > > MM,
	      unsigned long update_no,
	      XMLWriter& xml_out);

  private:
    InlineEigBndsMdagMParams params;
  };

};

#endif
