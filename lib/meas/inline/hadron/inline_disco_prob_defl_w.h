// -*- C++ -*-
/*! \file
 * \brief Compute disconnected diagrams with 4d hadamard probing and deflation
 *   3D probing is also allowed (I think ...)
 *
 * Propagator calculation on a colorvector
 */

#ifndef __inline_hada_disco_4d_w_h__
#define __inline_hada_disco_4d_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineDiscoProbDefl 
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
	int max_path_length ; /*! maximum displacement path */
	int p2_max ; /*! maximum p2  */
        multi2d<int> p_list; //Instead of a max momentum, a list is possible as an input.
        int p_num;  //Maximum number of momenta in the file.
	std::string p_file; //Name of file that contains list of momenta.
        bool use_p_list; //A boolean that keeps track of which momentum structure to pass to the fourier transform.
	bool multifile_write; //A boolean that switches between new and old code for writing to multiple databases,
	std::string mass_label ; /*! a std::string flag maybe used in analysis*/
        int max_rhs;            /*! maximum number of linear systems solved simultaneously */
	ChromaProp_t prop;
        GroupXML_t projParam;
	int probing_distance;
        std::string probing_file;
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
