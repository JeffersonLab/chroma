// $Id: qprop_io.h,v 1.8 2004-04-15 14:43:24 bjoo Exp $
/*! \file
 * \brief Routines associated with Chroma propagator IO
 */

#ifndef __qprop_io_h__
#define __qprop_io_h__

#include "handle.h"
#include "io/param_io.h"
#include "io/fermact_paramio.h"
#include "meas/sources/srcsnktype.h"
#include "meas/sources/wavetype.h"
#include "meas/smear/wvfkind.h"

/*
 * Chroma propagator support
 *
 * \ingroup io
 *
 * @{
 */


//! Propagator source header
struct PropSource_t
{
  int              version;
  SourceType       source_type;   // Point, Shell, Wall, etc.
  WaveStateType    wave_state;    // S-wave or P-wave
  SmearingParam_t  sourceSmearParam;
                             // wvf-function smearing type (Gaussian, Exponential, etc.)
                             // smearing width
                             // number of iteration for smearing
  int              j_decay;         // Decay direction
  int              direction;       // S-wave;   P-wave x(0) y(1) z(2)
  int              laplace_power;   // power=1 implies 1 laplacian
  int              disp_length;     // displacement length
  int              disp_dir;        // x(0), y(1), z(2)
  Real             link_smear_fact; // smearing factor
  int              link_smear_num;  // number of smearing hits
  multi1d<int>     t_source;        // source location
  multi1d<int>     nrow;            // lattice size
};


//! Propagator sink header
struct PropSink_t
{
  int              version;
  SinkType         sink_type;       // Point, Shell, Wall, etc.
  WaveStateType    wave_state;      // S-wave or P-wave
  SmearingParam_t  sinkSmearParam;
                             // wvf-function smearing type (Gaussian, Exponential, etc.)
                             // smearing width
                             // number of iteration for smearing
  int              direction;       // S-wave;   P-wave x(0) y(1) z(2)
  int              laplace_power;   // power=1 implies 1 laplacian
  int              disp_length;     // displacement length
  int              disp_dir;        // x(0), y(1), z(2)
  Real             link_smear_fact; // smearing factor
  int              link_smear_num;  // number of smearing hits
  multi1d<int>     nrow;            // lattice size
};


//! Propagator inversion parameters
class ChromaProp_t { 
 public:
  int             version;
  FermType        FermTypeP;
  
  // Uninitialised to start with (should hold null pointer)
  FermActParams* FermActHandle;
  
  InvertParam_t   invParam;   // Inverter parameters
  multi1d<int>    boundary;
  multi1d<int>    nrow;          // lattice size
  
  ChromaProp_t(void) { FermActHandle = 0x0 ; }
  ~ChromaProp_t(void) { if ( FermActHandle != 0x0 ) delete FermActHandle ; }
};


//! Structure for writing to seqprop files
struct ChromaSeqProp_t
{
  bool             nonRelSeqProp;
  InvertParam_t    invParam;
  int              Seq_src;
  multi1d<int>     sink_mom;
  int              t_sink;
  multi1d<int>     nrow;
};




//! Initialize header with default values
void initHeader(PropSource_t& header);

//! Initialize header with default values
void initHeader(PropSink_t& header);

//! Initialize header with a source header
void initHeader(PropSink_t& header, const PropSource_t& source);

//! Initialize header with default values
void initHeader(ChromaProp_t& header);

//! Propagator source read
void read(XMLReader& xml, const std::string& path, PropSource_t& header);

//! Propagator source writer
void write(XMLWriter& xml, const std::string& path, const PropSource_t& header);


//! Propagator sink reader
void read(XMLReader& xml, const std::string& path, PropSink_t& header);

//! Propagator sink writer
void write(XMLWriter& xml, const std::string& path, const PropSink_t& header);


//! Propagator header read
void read(XMLReader& xml, const std::string& path, ChromaProp_t& header);

//! Propagator header writer
void write(XMLWriter& xml, const std::string& path, const ChromaProp_t& header);


//! SeqPropagator header read
void read(XMLReader& xml, const std::string& path, ChromaSeqProp_t& header);

//! SeqPropagator header writer
void write(XMLWriter& xml, const std::string& path, const ChromaSeqProp_t& header);


//! Write a Chroma propagator
/*
 * \param file_xml     file header ( Read )
 * \param record_xml   xml holding propagator info ( Read )
 * \param quark_prop   propagator ( Read )
 * \param file         path ( Read )
 * \param volfmt       either QDP_SINGLEFILE, QDP_MULTIFILE ( Read )
 * \param serpar       either QDP_SERIAL, QDP_PARALLEL ( Read )
 */    
void writeQprop(XMLBufferWriter& file_xml,
		XMLBufferWriter& record_xml, const LatticePropagator& quark_prop,
		const string& file, 
		QDP_volfmt_t volfmt, QDP_serialparallel_t serpar);


//! Read a Chroma propagator
/*
 * \param file_xml     file header ( Write )
 * \param record_xml   xml holding propagator info ( Write )
 * \param quark_prop   propagator ( Write )
 * \param file         path ( Read )
 * \param serpar       either QDP_SERIAL, QDP_PARALLEL ( Read )
 */    
void readQprop(XMLReader& file_xml,
	       XMLReader& record_xml, LatticePropagator& quark_prop,
	       const string& file, 
	       QDP_serialparallel_t serpar);

/*! @} */  // end of group io

#endif
