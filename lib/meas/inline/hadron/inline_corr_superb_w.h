// -*- C++ -*-
/*! \file
 * \brief Inline measurement of baryon
 operators via colorstd::vector matrix elements
 */

#ifndef __inline_corr_superb_h__
#define __inline_corr_superb_h__

#include "io/xml_group_reader.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include <list>

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineCorrSuperbEnv 
  {
    bool registerAll();

    // Momentum list
    using mom_list_t = std::vector<std::vector<int>>;

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
	// Flavor to mass entry
	struct FlavorToMass {
	  char flavor;
	  std::string mass;
	};

	// Flavor to mass entry
	struct FlavorToProp {
	  char flavor;
	  ChromaProp_t prop;
	};

	int 			num_vecs;               /*! rank of the distillation basis */
	int 			max_rhs;                /*! maximum rhs to solve at once */
	int                     decay_dir;              /*!< Decay direction */
	std::vector<FlavorToMass> flavor_to_mass;       /*!< map from flavor to mass label */
	std::vector<FlavorToProp> flavor_to_prop;       /*!< map from flavor to prop */
	int			t_origin;		/*!< t_origin */
	int			Nt_forward;		/*!< Nt_forward */
	GroupXML_t              link_smearing;          /*!< link smearing xml */
	std::vector<std::string>  meson_files;          /*!< list of meson files */
	std::vector<std::string>  baryon_files;         /*!< list of baryons files */
	std::vector<std::string>  prop_files;           /*!< list of props files */
	std::vector<std::string>  genprop_files;        /*!< list of genprops files */
	std::string             ensemble;               /*!< ensemble name */
	bool                      testing;              /*!< whether to do tests */
      };

      struct NamedObject_t
      {
	std::string                 gauge_id;           /*!< Gauge field */
	std::vector<std::string>    colorvec_files;     /*!< Eigenvectors in mod format */
	std::string                 corr_graph_file;    /*!< Input corr graph */
	std::string                 corr_file;          /*!< File name where to save the corr functions */
      };

      Param_t        param;      /*!< Parameters */    
      NamedObject_t  named_obj;  /*!< Named objects */
      std::string    xml_file;   /*!< Alternate XML file pattern */
    };


    //! Inline measurement of baryon operators via colorstd::vector matrix elements
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

  } // namespace InlineCorrSuperbEnv 
}

#endif
