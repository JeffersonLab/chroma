// -*- C++ -*-
// $Id: inline_ritz_H_w.h,v 3.2 2007-08-23 19:02:44 edwards Exp $
/*! \file
 * \brief Inline construction of eigenvalues (Ritz)
 *
 * Eigenvalue calculations
 */

#ifndef __inline_ritz_H_w_h__
#define __inline_ritz_H_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/eigen_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineRitzEnv 
  {
    extern const std::string name;
    bool registerAll();


    //! Parameter structure
    /*! \ingroup inlinehadron */ 
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path);

      unsigned long     frequency;

      struct Param_t
      {
	int             version;
	GroupXML_t      fermact;          /*!< fermion action */
	RitzParams_t    ritz_params;
      } param;
      std::string       stateInfo;
    
      struct NamedObject_t
      {
	std::string     gauge_id;
	std::string     eigen_id;
      } named_obj;

      std::string xml_file;  // Alternate XML file pattern
    };

    //! Inline measurement of eigenvalues
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

  } // namespace InlineRitzEnv

} // namespace Chroma
#endif
