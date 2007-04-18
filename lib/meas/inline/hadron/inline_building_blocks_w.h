// -*- C++ -*-
// $Id: inline_building_blocks_w.h,v 3.5 2007-04-18 02:32:26 edwards Exp $
/*! \file
 * \brief Inline construction of BuildingBlocks
 *
 * Building Blocks on forward and sequential props
 */

#ifndef __inline_building_blocks_h__
#define __inline_building_blocks_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineBuildingBlocksEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineBuildingBlocksParams 
  {
    InlineBuildingBlocksParams();
    InlineBuildingBlocksParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    //! Parameters
    struct Param_t
    {
      bool         use_sink_offset;    /*!< should insertion origin be sink_mom */
      int          mom2_max;           /*!< (mom)^2 <= mom2_max */
      int          links_max;          /*!< maximum number of links */
      bool         canonical;          /*!< True if mom in BB filenames is canonicalized */
      bool         time_reverse;       /*!< Time reverse the building blocks */
      bool         translate;          /*!< Shifts the BB correlator output to start at t_source as 0*/
      GroupXML_t   cfs;                /*!< Fermion state */
    } param;

    //! Propagators
    struct NamedObject_t
    {
      std::string   BkwdPropId;         /*!< backward propagator */
      std::string   BkwdPropG5Format;   /*!< backward propagators Gamma5 Format */
      int           GammaInsertion;     /*!< second gamma insertion */
      std::string   Flavor;             /*!< Flavor id - values like U, D, S, C, B, T */
      std::string   BBFileNamePattern;  /*!< BB output file name pattern */
    };

    //! BB output
    struct BB_out_t
    {
      std::string       OutFileName;
      std::string       GaugeId;        /*!< Input Gauge id */
      std::string       FrwdPropId;     /*!< Input forward prop */
      multi1d<NamedObject_t>   BkwdProps;
    } bb;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of 3-pt functions writing building-blocks
  /*! \ingroup inlinehadron */
  class InlineBuildingBlocks : public AbsInlineMeasurement 
  {
  public:
    ~InlineBuildingBlocks() {}
    InlineBuildingBlocks(const InlineBuildingBlocksParams& p) : params(p) {}
    InlineBuildingBlocks(const InlineBuildingBlocks& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineBuildingBlocksParams params;
  };

};

#endif
