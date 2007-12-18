// -*- C++ -*-
// $Id: dilution_quark_source_const_w.h,v 1.3 2007-12-18 13:40:25 edwards Exp $
/*! \file
 * \brief Dilution operator constructed from pre-generated solutions.
 * 
 * This dilution operator takes pre-generated LatticeFermion solutions built
 * in the chain of MAKE_SOURCE and PROPAGATOR inline measurement calls
 */

#ifndef __dilution_quark_source_const_h__
#define __dilution_quark_source_const_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "meas/hadron/dilution_operator.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  /*! \ingroup inlinehadron */
  namespace DilutionQuarkSourceConstEnv 
  {
    extern const std::string name;
    bool registerAll();

    //! Parameter structure
    /*! \ingroup inlinehadron */
    struct Params
    {
      Params();
      Params(XMLReader& xml_in, const std::string& path);
      void writeXML(XMLWriter& xml_out, const std::string& path) const;

      //! Solution files for each quark
      struct QuarkFiles_t
      {
	//! Time dilution components  
	struct TimeDilutions_t
	{
	  multi1d<std::string> dilution_files;  /*!< dilution files for this time dilution*/
	};
	
	multi1d<TimeDilutions_t> time_files;
      };

      QuarkFiles_t  quark;       /*!< All the solutions that are needed for a single quark */
    }; // struct Params


    //! Structure holding a source and its solutions
    struct QuarkSourceSolutions_t
    {
      //! Structure holding solutions
      struct TimeDilutions_t
      {
	struct Dilutions_t
	{
	  LatticeFermion     source;
	  LatticeFermion     soln;
	  PropSourceConst_t  source_header;
	  ChromaProp_t       prop_header;
	};
	
	multi1d<int>       t0;  //Times included in this dilution
				
	multi1d<Dilutions_t>  dilutions; //dilutions per timeslice per spin
      };
				
      int   decay_dir;
      Seed  seed;
      multi1d<TimeDilutions_t>  time_dilutions;
    };


    //! Dilutions constructed from propagator solutions over MAKE_SOURCE dilutions
    class Dilute : public DilutionOperator<LatticeFermion>
    {
    public:
      //! Construct from a parameters struct
      Dilute(const Params& p) : params(p) {init();}

      //! Virtual destructor to help with cleanup;
      ~Dilute() {}

      //! Get an iterator for this dilution
      const_iterator begin() const;

      //! Get an iterator for this dilution
      const_iterator end() const;

      //! The decay direction
      int getDecayDir() const {return quark.decay_dir;}

      //! The seed identifies this quark
      const Seed& getSeed() const {return quark.seed;}

      //! Does this creation operator have support on times slice t0
      bool hasTimeSupport(const_iterator iter, int t0) const;

      //! Return the original source vector
      /*! NOTE: this might be slow */
      LatticeFermion source() const;
    
      //! Return the diluted source vector
      LatticeFermion dilutedSource(const_iterator iter) const;
    
      //! Return the solution vector corresponding to the diluted source
      LatticeFermion dilutedSolution(const_iterator iter) const;

    protected:
      //! Initialize the object
      void init();

    private:
      //! Hide default constructor
      Dilute() {}

    private:
      Params params;
      QuarkSourceSolutions_t quark;
    };
    
  } // namespace DilutionQuarkSourceConstEnv


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, DilutionQuarkSourceConstEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const DilutionQuarkSourceConstEnv::Params& param);

} // namespace Chroma

#endif
