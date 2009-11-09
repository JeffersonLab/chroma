// -*- C++ -*-
// $Id: inline_eigbnds.h,v 3.4 2008-01-30 18:29:07 bjoo Exp $

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
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineeig */
  namespace InlineEigBndsMdagMEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Parameter structure
  /*! \ingroup inlineeig */
  struct InlineEigBndsMdagMParams 
  {
    InlineEigBndsMdagMParams();
    InlineEigBndsMdagMParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_in, const std::string& path);

    unsigned long   frequency;
    GroupXML_t      ferm_act;    /*!< fermion action */
    bool            usePV;       /*!< measure eigs of PV matrix if applicable */

    //! Struct for parameters needed for a Ritz type solve
    struct RitzParams_t
    {
      Real RsdR;
      Real RsdA;
      Real RsdRHi;
      Real RsdAHi;
      Real RsdZero;
      bool ProjApsiP;
      int  Nmin;
      int  MaxCG;
      int  Nrenorm;
      int Neig;
    } ritz;

    struct NamedObject_t
    {
      std::string   gauge_id;
    } named_obj; 

    std::string xml_file;  // Alternate XML file pattern
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

    //! Do the measurement
    void operator()(unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

    //! 4D
    void do4d(Handle< LinearOperator<LatticeFermion> > MM,
	      unsigned long update_no,
	      XMLWriter& xml_out);

    //! 5D
    void do5d(Handle< LinearOperatorArray<LatticeFermion> > MM,
	      unsigned long update_no,
	      XMLWriter& xml_out);

  private:
    InlineEigBndsMdagMParams params;
    Handle< FermionAction<LatticeFermion,
			  multi1d<LatticeColorMatrix>,
			  multi1d<LatticeColorMatrix> > > fermact;
  };

}

#endif
