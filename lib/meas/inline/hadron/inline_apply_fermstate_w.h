// -*- C++ -*-
// $Id: inline_apply_fermstate_w.h,v 3.2 2006-11-21 04:07:46 edwards Exp $
/*! \file
 *  \brief Inline ferm state application
 */

#ifndef __inline_apply_fermstate_w_h__
#define __inline_apply_fermstate_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineFermStateEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  /*! \ingroup inlinehadron */
  struct InlineFermStateParams 
  {
    InlineFermStateParams();
    InlineFermStateParams(XMLReader& xml_in, const std::string& path);
    void writeXML(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      GroupXML_t    cfs;      /*!< Ferm State */
    } param;

    struct NamedObject_t
    {
      std::string   gauge_id;
      std::string   output_id;
    } named_obj;
  };


  /*! \ingroup inlinehadron */
  class InlineFermState : public AbsInlineMeasurement 
  {
  public:
    ~InlineFermState() {}
    InlineFermState(const InlineFermStateParams& p) : params(p) {}
    InlineFermState(const InlineFermState& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    void operator()(unsigned long update_no,
		    XMLWriter& xml_out); 

  private:
    InlineFermStateParams params;
  };

};

#endif
