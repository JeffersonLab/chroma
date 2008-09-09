// -*- C++ -*-
// $Id: inline_qqqNucNuc_w.h,v 3.5 2008-09-09 20:30:42 kostas Exp $
/*! \file
 * \brief The QQQ and QQBAR object calculation
 *
 */

#ifndef __inline_QQQNucNuc_h__
#define __inline_QQQNucNuc_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineQQQNucNucEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Parameter structure
  /*! \ingroup inlinehadron */
  struct InlineQQQNucNucParams 
  {
    InlineQQQNucNucParams();
    InlineQQQNucNucParams(XMLReader& xml_in, const std::string& path);
    void write(XMLWriter& xml_out, const std::string& path);

    unsigned long frequency;

    struct Param_t
    {
      int max_p2;
      bool doVectorMesons ;
      bool doDecupletBar  ;
    } param;

    struct NamedObject_t
    {
      std::string          gauge_id;  /*!< Input gauge field */
      multi1d<std::string> prop_ids;  /*!< Input forward propagators */
    } named_obj;

    std::string qqq_file ;  // binary file to write the qqq object
    std::string qqbar_file ;  // binary file to write the qqbar object
    
    std::string xml_file ;  // Alternate XML file pattern
  };


  //! Inline measurement of baryon-baryon 2-pt correlators
  /*! \ingroup inlinehadron */
  class InlineQQQNucNuc : public AbsInlineMeasurement 
  {
  public:
    ~InlineQQQNucNuc() {}
    InlineQQQNucNuc(const InlineQQQNucNucParams& p) : params(p) {}
    InlineQQQNucNuc(const InlineQQQNucNuc& p) : params(p.params) {}

    unsigned long getFrequency(void) const {return params.frequency;}

    //! Do the measurement
    void operator()(const unsigned long update_no,
		    XMLWriter& xml_out); 

  protected:
    //! Do the measurement
    void func(const unsigned long update_no,
	      XMLWriter& xml_out); 

  private:
    InlineQQQNucNucParams params;
  };

};

#endif
