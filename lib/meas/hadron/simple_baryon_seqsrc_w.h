// -*- C++ -*-
// $Id: simple_baryon_seqsrc_w.h,v 3.3 2006-11-28 20:00:49 edwards Exp $
/*! \file
 *  \brief Construct baryon sequential sources.
 */

#ifndef __simple_baryon_seqsrc_w_h__
#define __simple_baryon_seqsrc_w_h__

#include "meas/hadron/baryon_seqsrc_w.h"

namespace Chroma 
{

  //! Name and registration
  /*! @ingroup hadron */
  namespace SimpleBaryonSeqSourceEnv
  {
    bool registerAll();

  
    //! Simple baryon sequential source parameters
    /*! @ingroup hadron */
    struct Params
    {
      Params();
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;

      multi1d<int>     sink_mom;        /*!< sink momentum */
      int              t_sink;          /*!< time slice of sink */
      int              j_decay;         /*!< Decay direction */
    };


    //! Nucleon-Nucleon U piece with general projector and Cg5
    /*! @ingroup hadron
     *
     * Create a simple baryon sequential propagator source
     */
    class BarNuclUTCg5 : public BaryonSeqSourceBase
    {
    public:
      //! Full constructor
      BarNuclUTCg5(const Params& p, const SpinMatrix& spinT, const SpinMatrix& spinCg5) :
	params(p), T(spinT), Cg5(spinCg5) {}

      //! Default destructor
      ~BarNuclUTCg5() {}
      
      //! Construct the source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);

    protected: 
      //! Set bc
      multi1d<int>& getBC() {return bc;}

      //! Get bc
      const multi1d<int>& getBC() const {return bc;}

      //! Set t_srce
      multi1d<int>& getTSrce() {return t_srce;}

      //! Get t_srce
      const multi1d<int>& getTSrce() const {return t_srce;}

      //! Get t_sink
      int getTSink() const {return params.t_sink;}

      //! Get sink_mom
      const multi1d<int>& getSinkMom() const {return params.sink_mom;}

      //! Get decay_dir
      const int getDecayDir() const {return params.j_decay;}

    private:
      //! Hide partial constructor
      BarNuclUTCg5() {}

    private:
      multi1d<int>  t_srce;   /*<! Must come from propagator headers */
      multi1d<int>  bc;       /*<! Must come from propagator headers */
      Params  params;   /*!< Seqsource params */
      SpinMatrix T;     /*!< The spin projector matrix */
      SpinMatrix Cg5;   /*!< The Cg5 at the source and sink */
    };


    //! Nucleon-Nucleon D piece with general projector and Cg5
    /*! @ingroup hadron
     *
     * Create a simple baryon sequential propagator source
     */
    class BarNuclDTCg5 : public BaryonSeqSourceBase
    {
    public:
      //! Full constructor
      BarNuclDTCg5(const Params& p, const SpinMatrix& spinT, const SpinMatrix& spinCg5) :
	params(p), T(spinT), Cg5(spinCg5) {}

      //! Default destructor
      ~BarNuclDTCg5() {}
      
      //! Construct the source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);

    protected: 
      //! Set bc
      multi1d<int>& getBC() {return bc;}

      //! Get bc
      const multi1d<int>& getBC() const {return bc;}

      //! Set t_srce
      multi1d<int>& getTSrce() {return t_srce;}

      //! Get t_srce
      const multi1d<int>& getTSrce() const {return t_srce;}

      //! Get t_sink
      int getTSink() const {return params.t_sink;}

      //! Get sink_mom
      const multi1d<int>& getSinkMom() const {return params.sink_mom;}

      //! Get decay_dir
      const int getDecayDir() const {return params.j_decay;}

    private:
      //! Hide partial constructor
      BarNuclDTCg5() {}

    private:
      multi1d<int>  t_srce;   /*<! Must come from propagator headers */
      multi1d<int>  bc;       /*<! Must come from propagator headers */
      Params  params;   /*!< Seqsource params */
      SpinMatrix T;     /*!< The spin projector matrix */
      SpinMatrix Cg5;   /*!< The Cg5 at the source and sink */
    };


    //! Delta+ - Delta+ U piece with general projector and spin matrix
    /*! @ingroup hadron
     *
     * Create a simple baryon sequential propagator source
     */
    class BarDeltaUTsp : public BaryonSeqSourceBase
    {
    public:
      //! Full constructor
      BarDeltaUTsp(const Params& p, const SpinMatrix& spinT, const SpinMatrix& spinSP) :
	params(p), T(spinT), sp(spinSP) {}

      //! Default destructor
      ~BarDeltaUTsp() {}
      
      //! Construct the source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);

    protected: 
      //! Set bc
      multi1d<int>& getBC() {return bc;}

      //! Get bc
      const multi1d<int>& getBC() const {return bc;}

      //! Set t_srce
      multi1d<int>& getTSrce() {return t_srce;}

      //! Get t_srce
      const multi1d<int>& getTSrce() const {return t_srce;}

      //! Get t_sink
      int getTSink() const {return params.t_sink;}

      //! Get sink_mom
      const multi1d<int>& getSinkMom() const {return params.sink_mom;}

      //! Get decay_dir
      const int getDecayDir() const {return params.j_decay;}

    private:
      //! Hide partial constructor
      BarDeltaUTsp() {}

    private:
      multi1d<int>  t_srce;   /*<! Must come from propagator headers */
      multi1d<int>  bc;       /*<! Must come from propagator headers */
      Params  params;   /*!< Seqsource params */
      SpinMatrix T;     /*!< The spin projector matrix */
      SpinMatrix sp;    /*!< The spin at the source and sink */
    };

    //! Delta+ - Delta+ D piece with general projector and spin matrix
    /*! @ingroup hadron
     *
     * Create a simple baryon sequential propagator source
     */
    class BarDeltaDTsp : public BaryonSeqSourceBase
    {
    public:
      //! Full constructor
      BarDeltaDTsp(const Params& p, const SpinMatrix& spinT, const SpinMatrix& spinSP) :
	params(p), T(spinT), sp(spinSP) {}

      //! Default destructor
      ~BarDeltaDTsp() {}
      
      //! Construct the source
      LatticePropagator operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& forward_props);

      //! Compute the 2-pt at the sink
      Complex twoPtSink(const multi1d<LatticeColorMatrix>& u,
			const multi1d<ForwardProp_t>& forward_headers,
			const multi1d<LatticePropagator>& forward_props,
			int gamma_insertion);

    protected: 
      //! Set bc
      multi1d<int>& getBC() {return bc;}

      //! Get bc
      const multi1d<int>& getBC() const {return bc;}

      //! Set t_srce
      multi1d<int>& getTSrce() {return t_srce;}

      //! Get t_srce
      const multi1d<int>& getTSrce() const {return t_srce;}

      //! Get t_sink
      int getTSink() const {return params.t_sink;}

      //! Get sink_mom
      const multi1d<int>& getSinkMom() const {return params.sink_mom;}

      //! Get decay_dir
      const int getDecayDir() const {return params.j_decay;}

    private:
      //! Hide partial constructor
      BarDeltaDTsp() {}

    private:
      multi1d<int>  t_srce;   /*<! Must come from propagator headers */
      multi1d<int>  bc;       /*<! Must come from propagator headers */
      Params  params;   /*!< Seqsource params */
      SpinMatrix T;     /*!< The spin projector matrix */
      SpinMatrix sp;    /*!< The spin at the source and sink */
    };

  }  // end namespace


  //! Reader
  /*! @ingroup hadron */
  void read(XMLReader& xml, const string& path, SimpleBaryonSeqSourceEnv::Params& param);

  //! Writer
  /*! @ingroup hadron */
  void write(XMLWriter& xml, const string& path, const SimpleBaryonSeqSourceEnv::Params& param);

}  // end namespace Chroma


#endif
