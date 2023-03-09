// -*- C++ -*-
/*! \file
 * \brief Inline hadron spectrum calculations
 *
 * Hadron spectrum calculations
 */

#ifndef __inline_meson_da_h__
#define __inline_meson_da_h__

#include <map>
#include <vector>

#include "chromabase.h"

#include "meas/inline/abs_inline_measurement.h"
#include "util/ferm/meson_ops.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineMesDAEnv 
  {
    extern const std::string name;
    bool registerAll();

  }


  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineMesDAParams 
  {

    InlineMesDAParams();
    InlineMesDAParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;


    struct Param_t
    {
      bool time_rev;           // Use time reversal in baryon spectroscopy

      int mom2_max;            // (mom)^2 <= mom2_max. mom2_max=7 in szin.
      bool avg_equiv_mom;      // average over equivalent momenta
      bool pz_only;           //  restrict momentum in the z direction
      int  pz_max;           //  restrict momentum in the z direction

      multi1d<MesonOps::State_t> states ;        // holds the states

      int gamma ;
      int z_max ;
      int z_dir ;
      
      std::string ensemble ; // a std::string describing this ensemble 
    } param;

    struct NamedObject_t
    {
      std::string  gauge_id;           /*!< Input gauge field */

      struct Props_t
      {
	std::string  up_id;
	std::string  down_id;
	std::string  strange_id;
	std::string  charm_id;
	std::string  bottom_id;

	std::string  src;
	std::string  snk;
      };

      multi1d<Props_t> props;

    } named_obj;

    std::string xml_file;  // Alternate XML file pattern

  };


  //! Inline measurement of hadron spectrum
  /*! \ingroup inlinehadron */
  class InlineMesDA : public AbsInlineMeasurement 
  {
  public:
    ~InlineMesDA() {}
    InlineMesDA(const InlineMesDAParams& p) : params(p) {}
    InlineMesDA(const InlineMesDA& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineMesDAParams params;
  };

};

#endif
