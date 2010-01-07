// -*- C++ -*-
// $Id: inline_create_colorvecs.h,v 3.3 2009-04-11 03:32:46 edwards Exp $
/*! \file
 * \brief Construct colorvectors via power iteration of the laplacian
 */

#ifndef __inline_create_colorvecs_h__
#define __inline_create_colorvecs_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineCreateColorVecsEnv 
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
	int         num_iter;    /*!< Number of hits/iterations of the gaussian smearing */
	Real        width;       /*!< Smearing width - same conventions as gaussian quark smearing */
	int         num_orthog;  /*!< Number of hits/iterations of orthogonalization step */
	GroupXML_t  link_smear;  /*!< link smearing xml */
      };

      struct NamedObject_t
      {
	std::string     gauge_id;      /*!< Gauge field */
	std::string     colorvec_id;   /*!< Output colorvec id */
	GroupXML_t      colorvec_obj;  /*!< Output colorvecs */
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

  } // namespace CreateColorVecsEnv

}

#endif
