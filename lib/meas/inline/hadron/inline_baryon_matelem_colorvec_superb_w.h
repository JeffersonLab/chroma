// -*- C++ -*-
/*! \file
 * \brief Inline measurement of baryon
 operators via colorstd::vector matrix elements
 */

#ifndef __inline_baryon_matelem_colorvec_superb_h__
#define __inline_baryon_matelem_colorvec_superb_h__

#include "io/xml_group_reader.h"
#include "meas/inline/abs_inline_measurement.h"
#include <list>

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace InlineBaryonMatElemColorVecSuperbEnv 
  {
    bool registerAll();


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
	struct Displacement_t
	{
	  multi1d<int>          left;                   /*!< Left quark displacement array */
	  multi1d<int>          middle;                 /*!< Middle quark displacement array */
	  multi1d<int>          right;                  /*!< Rifht quark displacement array */
	};

	bool                    use_derivP;		/*!< Use derivatives */
	int                     mom2_max;               /*!< (mom)^2 <= mom2_max */
	int                     displacement_length;    /*!< Displacement length for creat. and annih. ops */
	int                     num_vecs;               /*!< Number of color vectors to use */
	int                     decay_dir;              /*!< Decay direction */
	multi1d<Displacement_t> displacement_list;      /*!< Array of displacements list to generate */
	GroupXML_t              link_smearing;          /*!< link smearing xml */
	int			Nt_forward;		/*!< Nt_forward */
	int			t_source;		/*!< t_source */
	multi1d<int> 		t_slices; 		/*!< alternative to Nt_forward and t_source */
	multi1d<float>          phase;         		/*!< Phase to apply to colorvecs */
	int 			max_tslices_in_contraction; /*! maximum number of contracted tslices simultaneously */
	int 			max_moms_in_contraction;/*! maximum number of contracted momenta simultaneously */
      };

      struct NamedObject_t
      {
	std::string                 gauge_id;           /*!< Gauge field */
	std::vector<std::string>    colorvec_files;     /*!< Eigenvectors in mod format */
	std::string                 baryon_op_file;     /*!< File name for creation operators */
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

  } // namespace InlineBaryonMatElemColorVecSuperbEnv 
}

#endif
