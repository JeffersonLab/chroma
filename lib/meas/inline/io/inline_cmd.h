// -*- C++ -*-
/*! \file
 * \brief Inline task to execute commands
 */

#ifndef __inline_cmd_h__
#define __inline_cmd_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineExecuteCmdEnv 
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
	  std::vector<std::string>           cmd;       /*!< command lines*/
	  bool           only_on_master;      /*!< whether to execute all commands on master node */
	  int            max_attempts;        /*!< maximum number of times to try to launch a command */
	  bool           constrain_to_gpu;    /*!< whether to execute the command on the current GPU */
      };

      Param_t           param;
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

      //! Do the measurement
      void func(const unsigned long update_no,
		XMLWriter& xml_out); 

    private:
      Params params;
    };
  }
}

#endif
