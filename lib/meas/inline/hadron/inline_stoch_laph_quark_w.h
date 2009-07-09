// -*- C++ -*-
// $Id: inline_stoch_laph_quark_w.h,v 3.1 2009-07-09 02:13:21 jbulava Exp $
/*! \file
 * \brief Compute the laph-diluted sources and solutions. Write them out to a single db file.  
 *
 * Propagator calculation on a laph diluted source 
 */

#ifndef __inline_stoch_laph_quark_w_h__
#define __inline_stoch_laph_quark_w_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineLaphDilutedPropsEnv 
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
				int             decay_dir;      /*!< Decay direction */

				bool spin_dil; 									/*!< include full spin dilution? */

				multi1d< multi1d<int> > eig_vec_dils; /*!< the array of eigenvector dilutions */ 
				
				multi1d<int>    t_sources;      /*!< Array of time slice sources for props, assumes full time dilution */
	
				
				multi1d<Seed> ran_seeds;        /*!<Array of random seeds, this will be the number of unique noise sources created. */
				
				std::string     mass_label;     /*!< Some kind of mass label for these quarks */

				ChromaProp_t    prop;


      } param;

      struct NamedObject_t
      {
				std::string     gauge_id;     /*!< Gauge field */
				std::string     eigvec_id;    /*!< LatticeColorVector EigenInfo */
				std::string     prop_file; /*!< File name for propagators and sources */
      
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
	

  } // namespace PropMatElemColorVec

}

#endif
