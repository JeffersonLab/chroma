// -*- C++ -*-
// $Id: inline_multipole_w.h,v 3.2 2007-04-18 02:32:26 edwards Exp $
/*! \file
 *  \brief Inline multipole measurements
 */

#ifndef __inline_multipole_h__
#define __inline_multipole_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineMultipoleEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineMultipoleParams 
  {
    InlineMultipoleParams();
    InlineMultipoleParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      int max_L;
    } param;

    //! Propagators
    struct NamedObject_t
    {
      std::string       seqprop_id;         // sequential prop
      int               GammaInsertion;     // second gamma insertion
    };

    //! Multipole output
    struct Multipole_out_t
    {
      std::string             gauge_id;     /*!< input forward prop */
      std::string             prop_id;      /*!< input forward prop */
      multi1d<NamedObject_t>  seqprops;
    } pole;

    std::string xml_file;  // Alternate XML file pattern
  };

  //! Inline measurement of multipole expansion
  /*! \ingroup inlinehadron */
  class InlineMultipole : public AbsInlineMeasurement 
  {
  public:
    ~InlineMultipole() {}
    InlineMultipole(const InlineMultipoleParams& p) : params(p) {}
    InlineMultipole(const InlineMultipole& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineMultipoleParams params;
  };

};

#endif
