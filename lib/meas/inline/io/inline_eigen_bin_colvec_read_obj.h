// -*- C++ -*-
// $Id: inline_eigen_bin_colvec_read_obj.h,v 3.1 2008-06-18 21:38:28 edwards Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object reading
 */

#ifndef __inline_eigen_bin_colvec_read_obj_h__
#define __inline_eigen_bin_colvec_read_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineEigenBinColVecReadNamedObjEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlineio */
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long frequency;

      struct NamedObject_t
      {
	std::string   object_id;
      } named_obj;

      struct File_t
      {
	std::string   file_name;
      } file;
    };

    //! Inline reading of latticecolorvectors that are eigenvectors
    /*! \ingroup inlineio */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const Params& p) : params(p) {}
      InlineMeas(const InlineMeas& p) : params(p.params) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the writing
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      Params params;
    };

  }
}

#endif
