// -*- C++ -*-
// $Id: simple_baryon_2pt_w.h,v 1.3 2007-11-30 06:04:23 kostas Exp $
/*! \file
 *  \brief Construct simple baryon 2pt correlators
 */

#ifndef __simple_baryon_2pt_w_h__
#define __simple_baryon_2pt_w_h__

#include "meas/hadron/baryon_2pt_w.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace SimpleBaryon2PtEnv
  {
    bool registerAll();

  
    //! Simple baryon 2pt parameters
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      multi1d<int>     sink_mom;        /*!< sink momentum */
      int              decay_dir;       /*!< Decay direction */
    };


    //! Nucleon-Nucleon 2pt with general projector and Cg5
    /*! @ingroup hadron
     *
     * Create a simple baryon 2pt correlator
     */
    class BarNuclTCg5 : public Baryon2PtBase
    {
    public:
      //! Full constructor
      BarNuclTCg5(const Params& p, const SpinMatrix& spinT, const SpinMatrix& spinCg5) :
	params(p), T(spinT), Cg5(spinCg5) {}

      //! Default destructor
      ~BarNuclTCg5() {}
      
      //! Construct the correlators
      /*! Default implementation supplied */
      multi1d<Hadron2PtContraction_t> operator()(const multi1d<LatticeColorMatrix>& u);

    protected: 
      //! Set bc
      multi1d<int>& getBC() {return bc;}

      //! Get bc
      const multi1d<int>& getBC() const {return bc;}

      //! Set t_srce
      multi1d<int>& getTSrce() {return t_srce;}

      //! Get t_srce
      const multi1d<int>& getTSrce() const {return t_srce;}

      //! Get decay_dir
      int& getDecayDir() {return decay_dir;}

      //! Get decay_dir
      int getDecayDir() const {return decay_dir;}

    private:
      //! Hide partial constructor
      BarNuclTCg5() {}

    private:
      multi1d<int>  t_srce;   /*<! Must come from propagator headers */
      multi1d<int>  bc;       /*<! Must come from propagator headers */
      Params  params;   /*!< Seqsource params */
      SpinMatrix T;     /*!< The spin projector matrix */
      SpinMatrix Cg5;   /*!< The Cg5 at the source and sink */
    };


    //! Delta+ - Delta+ 2pt piece with general projector and spin matrix
    /*! @ingroup hadron
     *
     * Create a simple baryon sequential propagator source
     */
    class BarDeltaTsp : public Baryon2PtBase
    {
    public:
      //! Full constructor
      BarDeltaTsp(const Params& p, const SpinMatrix& spinT, const SpinMatrix& spinSP) :
	params(p), T(spinT), sp(spinSP) {}

      //! Default destructor
      ~BarDeltaTsp() {}
      
      //! Construct the source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

    protected: 
      //! Set bc
      multi1d<int>& getBC() {return bc;}

      //! Get bc
      const multi1d<int>& getBC() const {return bc;}

      //! Set t_srce
      multi1d<int>& getTSrce() {return t_srce;}

      //! Get t_srce
      const multi1d<int>& getTSrce() const {return t_srce;}

      //! Get decay_dir
      const int getDecayDir() const {return params.j_decay;}

    private:
      //! Hide partial constructor
      BarDeltaTsp() {}

    private:
      multi1d<int>  t_srce;   /*<! Must come from propagator headers */
      multi1d<int>  bc;       /*<! Must come from propagator headers */
      Params  params;   /*!< Seqsource params */
      SpinMatrix T;     /*!< The spin projector matrix */
      SpinMatrix sp;    /*!< The spin at the source and sink */
    };


    //! Delta+ - Delta+ 2pt piece with general projector and spin matrix
    /*! @ingroup hadron
     *
     * Create a simple baryon propagator with different source
     *  and sink diquarks
     */
    class BarDeltaTspSRCspSNK : public Baryon2PtBase
    {
    public:
      //! Full constructor
      BarDeltaTspSRCspSNK(const Params& p, const SpinMatrix& spinT, const SpinMatrix& spinSRC, const SpinMatrix& spinSNK, ) :
	params(p), T(spinT), spSRC(spinSRC),spSNK(spinSNK) {}

      //! Default destructor
      ~BarDeltaTsp() {}
      
      //! Construct the source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

    protected: 
      //! Set bc
      multi1d<int>& getBC() {return bc;}

      //! Get bc
      const multi1d<int>& getBC() const {return bc;}

      //! Set t_srce
      multi1d<int>& getTSrce() {return t_srce;}

      //! Get t_srce
      const multi1d<int>& getTSrce() const {return t_srce;}

      //! Get decay_dir
      const int getDecayDir() const {return params.j_decay;}

    private:
      //! Hide partial constructor
      BarDeltaTspSRCspSNK() {}

    private:
      multi1d<int>  t_srce;   /*<! Must come from propagator headers */
      multi1d<int>  bc;       /*<! Must come from propagator headers */
      Params  params;   /*!< Seqsource params */
      SpinMatrix T;     /*!< The spin projector matrix */
      SpinMatrix sp_src;    /*!< The spin at the source  */
      SpinMatrix sp_snk;    /*!< The sink at the source  */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, SimpleBaryon2PtEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const SimpleBaryon2PtEnv::Params& param);

}  // end namespace Chroma


#endif
