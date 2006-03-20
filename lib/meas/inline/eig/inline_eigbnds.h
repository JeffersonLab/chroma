// -*- C++ -*-
// $Id: inline_eigbnds.h,v 2.1 2006-03-20 04:22:02 edwards Exp $

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
  /*! \ingroup inlineeig */
  namespace InlineEigBndsMdagMEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Parameter structure
  /*! \ingroup inlineeig */
  struct InlineEigBndsMdagMParams 
  {
    InlineEigBndsMdagMParams();
    InlineEigBndsMdagMParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_in, const std::string& path);

    unsigned long   frequency;
    std::string     ferm_act;
    bool            usePV;   // measure eigs of PV matrix if applicable

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

    struct NamedObject_t
    {
      std::string   gauge_id;
    } named_obj;
  };


  //! Inline measurement of eigenvalue bounds of M^dag*M
  /*! \ingroup inlineeig */
  class InlineEigBndsMdagM : public AbsInlineMeasurement 
  {
  public:
    ~InlineEigBndsMdagM() {}
    InlineEigBndsMdagM(const InlineEigBndsMdagMParams& p);
    InlineEigBndsMdagM(const InlineEigBndsMdagM& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    void operator()(unsigned long update_no,
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
    Handle< const FermionAction<LatticeFermion> > fermact;
  };

};

#endif
