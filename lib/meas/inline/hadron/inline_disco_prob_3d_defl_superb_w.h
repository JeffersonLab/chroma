// -*- C++ -*-
/*! \file
 * \brief Compute disconnected diagrams with 4d hadamard probing and deflation
 *   3D probing is also allowed (I think ...)
 *
 * Propagator calculation on a colorvector
 */

#ifndef __inline_disco_prob_3d_defl_superb_w_h__
#define __inline_disco_prob_3d_defl_superb_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineDiscoProb3dDeflSuperb 
  {
    bool registerAll();

    // Displacement list
    using displacement_list_t = std::vector<std::vector<int>>;

    // Momentum list
    using mom_list_t = std::vector<std::vector<int>>;

    // Displacements-momenta combo
    struct displacements_moms_combo_t {
      displacement_list_t displacement_list;
      mom_list_t mom_list;
    };

    // Displacement-momenta combos
    using displacements_moms_combos_t = std::vector<displacements_moms_combo_t>;

    //! Parameter structure
    /*! \ingroup inlinehadron */ 
    struct Params 
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);

      unsigned long     frequency;

      struct Param_t
      {
	int max_path_length ; /*! maximum displacement path in the z direction */
	displacement_list_t alt_displacements;   /*!< Alternative displacement paths */
	int mom2_min;               /*!< (mom)^2 >= mom2_min */
	int mom2_max;               /*!< (mom)^2 <= mom2_max */
	mom_list_t mom_list;	    /*!< Alternative array of momenta to generate */
	displacements_moms_combos_t combos_list; /*!< Alternative list of disps and momenta */
	std::vector<int>  t_sources;      /*!< Array of time slice sources */
	std::string mass_label ; /*! a std::string flag maybe used in analysis*/
        int max_rhs;            /*! maximum number of linear systems solved simultaneously */
	ChromaProp_t prop;
        GroupXML_t projParam;
	int probing_distance;
	int probing_power;
	int first_color;
	int num_colors;
	int noise_vectors;
	bool use_ferm_state_links ;
	bool use_superb_format;  /*!< Whether use the superb format for storing the elementals */
	int num_vecs;
	int                       max_tslices_in_contraction;  /*! maximum number of contracted tslices simultaneously */
        GroupXML_t link_smearing;
      } param;

      struct NamedObject_t
      {
	std::string     gauge_id;    /*!< Gauge field */
	std::string     sdb_file;    /*!< the db file to store loops */
	std::string     defl_sdb_file;    /*!< the db file to store the deterministic part of the loops */
	std::string     ip_sdb_file;    /*!< the db file to store the angles between the projector and distillation basis */
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

  } // namespace PropColorVec


}

#endif
