// -*- C++ -*-
// $Id: inline_grid_prop_w.h,v 3.2 2008-11-28 05:13:58 kostas Exp $
/*! \file
 * \brief Compute the propagator elements    M^-1 * multi1d<LatticeColorVector>
 *
 * Propagator calculation on a colorvector
 */

#ifndef __inline_grid_prop_w_h__
#define __inline_grid_prop_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineGridPropEnv 
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

	
	struct Sources_t
	{
	  multi1d<int>                     spatial_mask_size;
	  multi1d<multi1d<multi1d<int> > > spatial_masks ;
	  int decay_dir;            /*!< Decay direction */
	  multi1d<int> t_sources;   /*!< Array of time slice sources for props */

	  bool smear ;
	  GroupXML_t  smr; /*!< xml holding smearing params */
	  GroupXML_t  displace; /*!< xml holding displacement params */
	  GroupXML_t  link_smear;  /*!< link smearing xml */


	};

	ChromaProp_t    prop ;
	Sources_t       src  ;
      } param;
      
      GroupXML_t        map_obj_params; /*!< Parameters for MapObj factory */
      
      struct NamedObject_t
      {
	std::string     gauge_id;       /*!< Gauge field */
	std::string     prop_id;        /*!< Id for output propagator solutions */
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

    void read(XMLReader& xml, const string& path, 
	      InlineGridPropEnv::Params::Param_t::Sources_t& input) ;
  } // namespace GridProp


}

#endif
