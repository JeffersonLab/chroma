// -*- C++ -*-
// $Id: inline_block_colorvecs.h,v 3.1 2009-02-23 17:03:15 kostas Exp $
/*! \file
 * \brief Compute the matrix element of  LatticeColorVector*M^-1*LatticeColorVector
 *
 * Propagator calculation on a colorvector
 */

#ifndef __inline_block_colorvecs_h__
#define __inline_block_colorvecs_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineBlockColorVecsEnv 
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

	struct Sources_t {
	  int decay_dir;            /*!< Decay direction */
	  GroupXML_t  smr; /*!< xml holding smearing params */
	  GroupXML_t  link_smear;  /*!< link smearing xml */
	};

	Sources_t       src  ;
       	bool OrthoNormal ; // if set to true will first do a global orthonormilization. Default is false
	multi1d<int> block ;
      } param;

      struct NamedObject_t
      {
	std::string     gauge_id;      /*!< Gauge field */
	std::string     colorvec_id;   /*!< Id for color vectors */
      } named_obj;

      std::string xml_file;  // Alternate XML file pattern
    };


    //! Inline task for compute LatticeColorVector matrix elements of a propagator
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

  } // namespace BlockColorVecsEnv

}

#endif
