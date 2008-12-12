// -*- C++ -*-
// $Id: inline_meson_grid_matelem_w.h,v 3.2 2008-12-12 03:54:57 kostas Exp $
/*! \file
 * \brief Inline measurement of meson operators via colorvector matrix elements
 */

#ifndef __inline_meson_grid_matelem_h__
#define __inline_meson_grid_matelem_h__

#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineMesonGridMatElemEnv 
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
	 struct Grid_t
	 {
	   multi1d<int>                     spatial_mask_size;
	   multi1d<multi1d<multi1d<int> > > spatial_masks ;
	   int decay_dir;            /*!< Decay direction */
	 };

	int                     mom2_max              ; /*!< (mom)^2 <= mom2_max                           */
	int                     displacement_length   ; /*!< Displacement length for creat. and annih. ops */
	multi1d< multi1d<int> > displacement_list     ; /*!< Array of displacements list to generate       */
	GroupXML_t              link_smearing         ; /*!< link smearing xml                             */
	GroupXML_t              smearing              ; /*!< quark smearing xml                            */
	Grid_t                  grid                  ; /*!< grid descriptor                               */
	bool                    smear                 ;
      };

      struct NamedObject_t
      {
	std::string         gauge_id;               /*!< Gauge field */
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

  } // namespace InlineMesonGridMatElemEnv 
}

#endif
