// -*- C++ -*-
/*! \file
 * \brief Inline task to read an object
 */

#ifndef __inline_read_subset_vectors_h__
#define __inline_read_subset_vectors_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineReadSubsetVectorsEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlineio */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned int frequency;

      struct File {
	std::string   file_name;
      } file;

      struct NamedObject_t {
	std::string   object_id;
      } named_obj;
    };

    //! Inline writing of map objects
    /*! \ingroup inlineio */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the writing
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      Params params;
    };

  } // namespace InlineReadSubsetVectorsEnv
}

#endif
