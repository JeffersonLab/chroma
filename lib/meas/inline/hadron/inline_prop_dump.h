// -*- C++ -*-
/*! \file
 * \brief Compute the propagator from distillation
 *
 * Propagator calculation in distillation
 */

#ifndef __inline_prop_dump_h__
#define __inline_prop_dump_h__

#ifndef QDP_IS_QDPJIT_NO_NVPTX

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlinePropDumpEnv 
  {
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlinehadron */ 
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long     frequency;

      struct Param_t
      {
	int           t_source;      /*!< Array of time slice sources for props */
	std::string   mass_label;     /*!< Some kind of mass label */
	int           num_vecs;      /*!< In case of bad things happening in the solution vectors, do retries */
	int           Nt_forward;      /*!< In case of bad things happening in the solution vectors, do retries */
	int           decay_dir;
      };

      struct NamedObject_t
      {
	std::string                 prop_op_file;    /*!< File name for propagator matrix elements */
      };

      Param_t           param;
      NamedObject_t     named_obj;
      std::string       xml_file;       /*!< Alternate XML file pattern */
    };


    //! Inline task for the propagator from distillation
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

  } // namespace PropColorVec


}

#endif // QDPJIT

#endif // Chroma
