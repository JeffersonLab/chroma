// -*- C++ -*-
/*! \file
 * \brief Inline measurement of tetra-quark operators via colorvector matrix elements
 */

#ifndef __inline_tetraquark_matelem_colorvec_h__
#define __inline_tetraquark_matelem_colorvec_h__

#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineTetraQuarkMatElemColorVecEnv 
  {
    bool registerAll();


    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path) const;

      unsigned long      frequency;

      struct Param_t
      {
	int                     mom2_min;               /*!< (mom)^2 >= mom2_min */
	int                     mom2_max;               /*!< (mom)^2 <= mom2_max */
	multi1d< multi1d<int> > mom_list;               /*!< Array of momenta to generate */
	int                     num_vecs;               /*!< Number of color vectors to use */
	#ifdef BUILD_SB
        bool use_superb_format = true;
#else
        bool use_superb_format = false;
#endif
	int                     decay_dir;              /*!< Decay direction */
      };

      struct NamedObject_t
      {
	std::string         colorvec_file;              /*!< LatticeColorVector EigenInfo */
	std::string         tetraquark_op_file;         /*!< File name for creation operators */
      };

      Param_t        param;      /*!< Parameters */    
      NamedObject_t  named_obj;  /*!< Named objects */
      std::string    xml_file;   /*!< Alternate XML file pattern */
    };


    //! Inline measurement of meson operators via colorvector matrix elements
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

  } // namespace InlineTetraQuarkMatElemColorVecEnv 
}

#endif
