// -*- C++ -*-
// $Id: inline_block_prop_matelem_w.h,v 1.1 2009-01-30 05:31:41 kostas Exp $
/*! \file
 * \brief Compute the matrix element of  LatticeColorVector*M^-1*LatticeColorVector
 *
 * Propagator calculation on a colorvector
 */

#ifndef __inline_block_prop_matelem_w_h__
#define __inline_block_prop_matelem_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineBlockPropMatElemEnv 
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
	int             num_vecs   ;  /*!< Number of color vectors to use */
	int             decay_dir  ;  /*!< Decay direction */
	multi1d<int>    t_sources  ;  /*!< Array of time slice sources for props */
	multi1d<int>    block      ;  /*!< describes the block */
	std::string     mass_label ;  /*!< Some kind of mass label */
      } param;

      struct NamedObject_t
      {
	std::string     gauge_id;       /*!< Gauge field */
	std::string     colorvec_id;    /*!< LatticeColorVector EigenInfo */
	std::string     prop_id;        /*!< Id for input propagator solutions */
	std::string     prop_op_file;   /*!< File name for propagator matrix elements */
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

  } // namespace BlockPropMatElem

}

#endif
