// -*- C++ -*-
/*! \file
 * \brief Compute disconnected diagrams with 4d hadamard probing and deflation
 *   3D probing is also allowed (I think ...)
 *
 * Propagator calculation on a colorvector
 */

#ifndef __inline_disco_frequencysplitting_prob_multimass_superb_w_h__
#define __inline_disco_frequencysplitting_prob_multimass_superb_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineDiscoFreqSplitProbMMSuperb
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
	int max_path_length ; 		/*! maximum displacement path */
	int p2_max ; 			/*! maximum p2  */
        multi2d<int> p_list; 		//Instead of a max momentum, a list is possible as an input.
        int p_num;  			//Maximum number of momenta in the file.
	std::string p_file; 		//Name of file that contains list of momenta.
        bool use_p_list; 		//A boolean that keeps track of which momentum structure to pass to the fourier transform.
	bool multifile_write; 		//A boolean that switches between new and old code for writing to multiple databases,
	std::string mass_label ; 	/*! a std::string flag maybe used in analysis*/
        int max_rhs;            	/*! maximum number of linear systems solved simultaneously*/
	int max_rhs_interp;		/* maxmum number of linear systems solved simulataneously during variance estimation */
        bool use_interpolation;         //Flag to do interpolation of the variances or not
	bool use_mg; 			//Flag to use the multigrid cost function or not
	int num_shifts;
        multi1d<Real> shifts;   	//The shifts used for frequency splitting. If use_interpolation == true, they are your test shifts
					//If use_interpolation == false, these will be the shifts you use for frequency splitting. 
	multi1d<double> del_s;		//The discretization for the interpolated shifts, in matlab the interpolated shifts are:  [0:del_s:shifts(end)];
	multi1d<int> num_bcshifts;
	ChromaProp_t prop;  		//propagator
	int probing_distance;
	int probing_power;
        std::string probing_file;	//probing file for the inverse
	multi1d<int> noise_vectors;  	//The number of noise vectors for each level for the trace estimation
	int num_samples;		//The number of noise vectors to estimate the variance for interpolation. Same number used per level
	bool use_ferm_state_links ;
	bool use_hpe;			//Flag to use hpe or not. This will be set to false if your shifts are not large enough for hpe
	std::vector<std::string> hpe_probing_files;  //Probing files to be used for exact computation of HPE trace
	int hpe_power;			//The power to truncate the HPE expansion
	multi1d<int> disp_gamma; 	//The gamma and displacement combination to optimize. If not included, one will be automatically chosen 
					//by finding the maximum relative error of non vanishing traces. 
	bool display_level_stats;	//Bool to display the stats for each level
	bool do_trace_est;
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
