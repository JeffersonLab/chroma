// -*- C++ -*-
// $Id: inline_sfpcac_w.h,v 1.3 2007-08-23 19:02:45 edwards Exp $
/*! \file
 * \brief Inline Schroedinger functional measurements
 */

#ifndef __inline_sfpcac_w_h__
#define __inline_sfpcac_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineSFpcacEnv 
  {
    extern const std::string name;
    bool registerAll();


    //! Parameter structure
    /*! \ingroup inlinehadron */ 
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long     frequency;

      ChromaProp_t      param;

      struct SFpcac_t
      {
	int  decay_dir;        /*!< decay direction */
	bool ZVfactP;          /*!< Measure Z_V */
	bool ZAfactP;          /*!< Measure Z_A */
	int  x0;               /*!< Starting location of currents */
	int  y0;               /*!< Ending location of currents */
      } sfpcac;

      struct NamedObject_t
      {
	std::string     gauge_id;
      } named_obj;

      std::string xml_file;  // Alternate XML file pattern
    };


    //! Inline measurement of Wilson loops
    /*! \ingroup inlinehadron */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the measurement
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    protected:
      //! Do the measurement
      void func(const unsigned long update_no,
		XMLWriter& xml_out); 

    private:
      Params params;
    };

  } // namespace InlineSFpcacEnv
}  // namespace Chroma

#endif
