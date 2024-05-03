// -*- C++ -*-
/*! \file
 * \brief Compute the propagator from distillation
 *
 * Propagator calculation in distillation
 */

#ifndef __inline_eigenvalues_superb_w_h__
#define __inline_eigenvalues_superb_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"

#ifdef BUILD_SB

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineEigenvaluesSuperbEnv 
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
	struct Contract_t
	{
	  int           num_vecs;       /*!< Number of eigenvectors to compute */
	  double        tol;            /*!< Tolerance of the eigenvectors to compute */
	  std::string   mass_label;     /*!< mass label */
          int           max_rhs;        /*! maximum number of linear systems solved simultaneously */
	};

	ChromaProp_t    prop;
	std::string     eigensolver; /* eigensolver options */
	Contract_t      contract;
      };

      struct NamedObject_t
      {
	std::string                 gauge_id;        /*!< Gauge field */
	std::string                 eigs_file;    /*!< File name for storing the eigenvalues and eigenvectors */
      };

      Param_t           param;
      NamedObject_t     named_obj;
      std::string       xml_file;       /*!< Alternate XML file pattern */
    };


    //! Inline task for the propagator from distillation
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

  } // namespace PropColorVec


}

#endif // BUILD_SB

#endif
