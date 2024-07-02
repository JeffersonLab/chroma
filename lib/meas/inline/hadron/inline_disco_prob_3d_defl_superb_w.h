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
	std::vector< std::vector<int>> alt_displacements;   /*!< Alternative displacement paths */
	int mom2_min;               /*!< (mom)^2 >= mom2_min */
	int mom2_max;               /*!< (mom)^2 <= mom2_max */
	std::vector<std::vector<int>>  mom_list;        /*!< Alternative array of momenta to generate */
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
      } param;

      struct NamedObject_t
      {
	std::string     gauge_id;    /*!< Gauge field */
	std::string     sdb_file;    /*!< the db file to store loops */
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
