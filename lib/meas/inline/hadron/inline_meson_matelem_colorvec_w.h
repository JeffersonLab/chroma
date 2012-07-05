// -*- C++ -*-
// $Id: inline_meson_matelem_colorvec_w.h,v 1.3 2008-08-06 17:58:00 edwards Exp $
/*! \file
 * \brief Inline measurement of meson operators via colorvector matrix elements
 */

#ifndef __inline_meson_matelem_colorvec_h__
#define __inline_meson_matelem_colorvec_h__

#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineMesonMatElemColorVecEnv 
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
	int                     displacement_length;    /*!< Displacement length for creat. and annih. ops */
	int                     num_vecs;               /*!< Number of color vectors to use */
	int                     decay_dir;              /*!< Decay direction */
	multi1d< multi1d<int> > displacement_list;      /*!< Array of displacements list to generate */
	GroupXML_t              link_smearing;          /*!< link smearing xml */

	// This all may need some work
	bool                    orthog_basis;           /*!< Whether all the basis vectors are orthog */
      };

      struct NamedObject_t
      {
	std::string         gauge_id;               /*!< Gauge field */
	std::string         colorvec_id;            /*!< LatticeColorVector EigenInfo */
	std::string         meson_op_file;          /*!< File name for creation operators */
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

  } // namespace InlineMesonMatElemColorVecEnv 
}

#endif
