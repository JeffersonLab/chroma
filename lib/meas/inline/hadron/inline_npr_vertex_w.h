// -*- C++ -*-
// $Id: inline_npr_vertex_w.h,v 1.4 2007-11-28 04:05:37 kostas Exp $
/*! \file
 * \brief Inline construction of NPR vertices
 *
 * NPR vertices on  props
 */

#ifndef __inline_npr_vertex_h__
#define __inline_npr_vertex_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineNprVertexEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineNprVertexParams 
  {
    InlineNprVertexParams();
    InlineNprVertexParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    //! Parameters
    struct Param_t
    {
      int          links_max;          /*!< maximum number of links */
      std::string  file_name;          /*!< bb output file name pattern */
      GroupXML_t   cfs;                /*!< Fermion state */
    } param;

    //! Propagators
    struct NamedObject_t
    {
      std::string       gauge_id;        /*!< Input Gauge id */
      std::string       prop_id;         /*!< Input forward prop */
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of NPR vertices
  /*! \ingroup inlinehadron */
  class InlineNprVertex : public AbsInlineMeasurement 
  {
  public:
    ~InlineNprVertex() {}
    InlineNprVertex(const InlineNprVertexParams& p) : params(p) {}
    InlineNprVertex(const InlineNprVertex& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineNprVertexParams params;
  };

};

#endif
