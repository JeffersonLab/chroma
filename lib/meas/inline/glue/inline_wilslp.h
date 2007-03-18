// -*- C++ -*-
// $Id: inline_wilslp.h,v 3.7 2007-03-18 21:59:15 edwards Exp $
/*! \file
 *  \brief Inline Wilson loops
 */

#ifndef __inline_wilslp_h__
#define __inline_wilslp_h__

#include "chromabase.h"
#include "handle.h"
#include "create_state.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlineglue */
  namespace InlineWilsonLoopEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Parameter structure
  /*! \ingroup inlineglue */
  struct InlineWilsonLoopParams 
  {
    InlineWilsonLoopParams();
    InlineWilsonLoopParams(XMLReader& xml_in, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      int           kind;
      int           j_decay;
      int           t_dir;
      GroupXML_t    cgs;      /*!< Gauge State */
    } param;

    struct NamedObject_t
    {
      std::string   gauge_id;
    } named_obj;

    std::string xml_file;  // Alternate XML file pattern
  };


  //! Inline measurement of Wilson loops
  /*! \ingroup inlineglue */
  class InlineWilsonLoop : public AbsInlineMeasurement 
  {
  public:
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    ~InlineWilsonLoop() {}
    InlineWilsonLoop(const InlineWilsonLoopParams& p) : params(p) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineWilsonLoopParams params;
  };

};

#endif
