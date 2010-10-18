// -*- C++ -*-
/*! \file
 * \brief Setup the origin and noise factory for distillution
 *
 * Setup for distillution
 */

#ifndef __inline_distillution_noise_h__
#define __inline_distillution_noise_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineDistillutionNoiseEnv 
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
	int             decay_dir;      /*!< Decay direction */
	std::string     ensemble;       /*!< Ensem label */
	std::string     sequence;       /*!< Sequence label - the trajectory */
      };

      struct NamedObject_t
      {
	std::string     distillution_id;/*!< Holds state for forming noises, including the origin */
      };

      Param_t           param;
      NamedObject_t     named_obj;
      std::string xml_file;  // Alternate XML file pattern
    };


    //! Inline task for setting up distillution
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

  } // namespace PropMatElemColorVec

}

#endif
