// -*- C++ -*-
/*! \file
 * \brief Construct colorvectors via power iteration of the laplacian
 */

#ifndef __inline_create_colorvecs_superb_h__
#define __inline_create_colorvecs_superb_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineCreateColorVecsSuperbEnv 
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
	int         num_vecs;    /*!< Number of vectors */
	int         decay_dir;   /*!< Decay direction */
	int         t_start;     /*!< First t_slice to compute */
	int         Nt_forward;  /*!< Number of t_slices to compute */
	GroupXML_t  link_smear;  /*!< link smearing xml */
	bool        write_fingerprint;  /*!< whether to write a portion of the colorvec instead of the whole */
	multi1d<float>  phase;   /*!< Phase to apply to colorvecs */
      };

      struct NamedObject_t
      {
	std::string     gauge_id;                           /*!< Gauge field */
	std::vector<std::string>    colorvec_files;         /*!< Eigenvectors in mod format */
	std::string     colorvec_out;                       /*!< Output colorvec */
      };

      Param_t        param;      /*!< Parameters */
      NamedObject_t  named_obj;  /*!< Named objects */
      std::string    xml_file;   /*!< Alternate XML file pattern */
    };


    //! Inline task for creating colorvectors
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

  } // namespace CreateColorVecsSuperbEnv

}

#endif
