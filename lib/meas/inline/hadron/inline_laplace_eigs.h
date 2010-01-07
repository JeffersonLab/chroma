// -*- C++ -*-
// $Id: inline_laplace_eigs.h,v 1.2 2009-06-23 15:12:42 jbulava Exp $
/*! \file
 * \brief Use the Implicitly Restarted Lanczos method with a Tchebyshev 
 * polynomial preconditioner to solve for the lowest eigenvalues and 
 * eigenvectors of the gague-covariant Laplacian
 */

#ifndef __inline_laplace_eigs_h__
#define __inline_laplace_eigs_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineLaplaceEigsEnv 
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
	int         num_vecs;    /*!< Number of vectors */
	int         decay_dir;   /*!< Decay direction */
	int         max_iter;    /*!< Maximum number of Lanczos iterations */
	Real 		tol; 		 /*!< Allowed residual upon exit */	

	GroupXML_t  link_smear;  /*!< link smearing xml */
      };

      struct NamedObject_t
      {
	std::string     gauge_id;      /*!< Gauge field */
	std::string     colorvec_id;   /*!< Id for color vectors */
	GroupXML_t      colorvec_obj;  /*!< Output colorvecs */
      };

      Param_t        param;      /*!< Parameters */
      NamedObject_t  named_obj;  /*!< Named objects */
      std::string    xml_file;   /*!< Alternate XML file pattern */
    };


    //! Inline task for Laplcian eigenvectors 
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

  } // namespace LaplaceEigsEnv

}

#endif
