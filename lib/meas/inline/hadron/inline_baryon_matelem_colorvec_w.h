// -*- C++ -*-
// $Id: inline_baryon_matelem_colorvec_w.h,v 3.3 2009-09-14 21:06:21 edwards Exp $
/*! \file
 * \brief Inline measurement of baryon
 operators via colorvector matrix elements
 */

#ifndef __inline_baryon_matelem_colorvec_h__
#define __inline_baryon_matelem_colorvec_h__

#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineBaryonMatElemColorVecEnv 
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
	struct Displacement_t
	{
	  multi1d<int>          left;                   /*!< Left quark displacement array */
	  multi1d<int>          middle;                 /*!< Middle quark displacement array */
	  multi1d<int>          right;                  /*!< Rifht quark displacement array */
	};

	bool                    use_derivP;             /*!< Are these displacements or derivatives? */
	int                     mom2_max;               /*!< (mom)^2 <= mom2_max */
	int                     displacement_length;    /*!< Displacement length for creat. and annih. ops */
	int                     num_vecs;               /*!< Number of color vectors to use */
	int                     decay_dir;              /*!< Decay direction */
	multi1d<Displacement_t> displacement_list;      /*!< Array of displacements list to generate */
	GroupXML_t              link_smearing;          /*!< link smearing xml */

	// This all may need some work
	bool                    site_orthog_basis;      /*!< Whether all the basis vectors are site level orthog */
      };

      struct NamedObject_t
      {
	std::string         gauge_id;               /*!< Gauge field */
	std::string         colorvec_id;            /*!< LatticeColorVector EigenInfo */
	std::string         baryon_op_file;          /*!< File name for creation operators */
      };

      Param_t        param;      /*!< Parameters */    
      NamedObject_t  named_obj;  /*!< Named objects */
      std::string    xml_file;   /*!< Alternate XML file pattern */
    };


    //! Inline measurement of baryon operators via colorvector matrix elements
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

  } // namespace InlineBaryonMatElemColorVecEnv 
}

#endif
