// -*- C++ -*-
/*! \file
 * \brief Inline task to write and delete an object from a named buffer
 *
 * Named object writing
 */

#ifndef __inline_qio_write_erase_obj_h__
#define __inline_qio_write_erase_obj_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "meas/inline/io/inline_qio_write_obj.h"

namespace Chroma 
{ 
  /*! \ingroup inlineio */
  namespace InlineQIOWriteEraseNamedObjEnv 
  {
    bool registerAll();

    //! Inline writing of memory objects
    /*! \ingroup inlineio */
    class InlineMeas : public AbsInlineMeasurement 
    {
    public:
      ~InlineMeas() {}
      InlineMeas(const InlineQIOWriteNamedObjEnv::Params& p) : params(p) {}

      unsigned long getFrequency(void) const {return params.frequency;}

      //! Do the writing
      void operator()(const unsigned long update_no,
		      XMLWriter& xml_out); 

    private:
      InlineQIOWriteNamedObjEnv::Params params;
    };

  }

}

#endif
