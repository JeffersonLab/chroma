// -*- C++ -*-
// $Id: inline_noisy_building_blocks_w.h,v 2.1 2005-10-21 15:55:39 flemingg Exp $
/*! \file
 * \brief Inline construction of noisy BuildingBlocks
 *
 * Noisy Building Blocks on forward props and noisy sources
 */

#ifndef __inline_noisy_building_blocks_h__
#define __inline_noisy_building_blocks_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineNoisyBuildingBlocksEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineNoisyBuildingBlocksParams 
  {
    InlineNoisyBuildingBlocksParams();
    InlineNoisyBuildingBlocksParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    //! Parameters
    struct Param_t
    {
      bool     use_sink_offset;    /*!< should insertion origin be sink_mom */
      int      mom2_max;           /*!< (mom)^2 <= mom2_max */
      int      links_max;          /*!< maximum number of links */
      bool     canonical;          /*!< True if mom in BB filenames is canonicalized */
      multi1d<int> nrow;           /*!< lattice size */
    } param;

    //! Propagators
//  struct NamedObject_t
//  {
//    std::string   BkwdPropId;         // backward propagator
//    std::string   BkwdPropG5Format;   // backward propagators Gamma5 Format
//    int           GammaInsertion;     // second gamma insertion
//    std::string   Flavor;             // Flavor id - values like U, D, S, C, B, T
//    std::string   BBFileNamePattern;  // BB output file name pattern
//  };

    //! BB output
    struct BB_out_t
    {
      std::string       OutFileName;
      std::string       FrwdPropId;     // input forward prop
//    multi1d<NamedObject_t>   BkwdProps;
      std::string       NoisySrcId;     // input noisy source
      std::string   NoisySrcG5Format;   // noisy source Gamma5 Format
      int           GammaInsertion;     // second gamma insertion
      std::string   Flavor;             // Flavor id - values like U, D, S, C, B, T
      std::string   BBFileNamePattern;  // BB output file name pattern
    } bb;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlinehadron */
  class InlineNoisyBuildingBlocks : public AbsInlineMeasurement 
  {
  public:
    ~InlineNoisyBuildingBlocks() {}
    InlineNoisyBuildingBlocks(const InlineNoisyBuildingBlocksParams& p) : params(p) {}
    InlineNoisyBuildingBlocks(const InlineNoisyBuildingBlocks& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const multi1d<LatticeColorMatrix>& u,
		    XMLBufferWriter& gauge_xml,
		    const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const multi1d<LatticeColorMatrix>& u,
	      XMLBufferWriter& gauge_xml,
	      const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineNoisyBuildingBlocksParams params;
  };

};

#endif
