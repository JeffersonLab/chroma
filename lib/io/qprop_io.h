// -*- C++ -*-
// $Id: qprop_io.h,v 1.27 2005-03-07 02:54:15 edwards Exp $
/*! \file
 * \brief Routines associated with Chroma propagator IO
 */

#ifndef __qprop_io_h__
#define __qprop_io_h__

#include "io/smearing_io.h"
#include "io/enum_io/enum_io.h"
#include "meas/sources/srcsnktype.h"
#include "meas/sources/wavetype.h"
#include "meas/smear/wvfkind.h"

namespace Chroma {

/*!
 * Chroma propagator support
 *
 * \ingroup io
 *
 * @{
 */


//! Propagator source header
struct PropSource_t
{
  PropSource_t();                 // default constructor
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
  PropSink_t();                     // default constructor
  PropSink_t(const PropSource_t& source); // from a prop source

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


struct ChromaMultiProp_t 
{ 
  ChromaMultiProp_t();          // default constructor
  bool            nonRelProp;   // compute only the nonrelativistic portion of this prop?

  // String holding XML of the FermionAction section
  std::string     fermact;

  MultiInvertParam_t   invParam;   // Inverter parameters
  multi1d<int>    nrow;          // lattice size
 
  multi1d<Real>   MultiMasses;
};


//! Propagator inversion parameters
struct ChromaProp_t 
{ 
  ChromaProp_t();               // default constructor
  bool            nonRelProp;   // compute only the nonrelativistic portion of this prop?

  // String holding XML of the FermionAction section
  std::string     fermact;
  
  // String holding XML for auxiliary state information
  InvertParam_t   invParam;   // Inverter parameters
  multi1d<int>    nrow;          // lattice size
};


//! Structure for writing to seqsource files
struct SeqSource_t
{
  std::string      seq_src;
  multi1d<int>     sink_mom;
  int              t_sink;
  multi1d<int>     nrow;
};


//! Structure for writing to seqprop files
struct ChromaSeqProp_t
{
  bool             nonRelSeqProp;
  InvertParam_t    invParam;
  std::string      seq_src;
  multi1d<int>     sink_mom;
  int              t_sink;
  multi1d<int>     nrow;
};


//! Mega structure holding a full forward prop (except gauge)
struct ForwardProp_t
{
  PropSink_t       sink_header;
  ChromaProp_t     prop_header;
  PropSource_t     source_header;
};


//! Mega structure holding a full sequential source (except gauge)
struct SequentialSource_t
{
  PropSink_t       sink_header;
  SeqSource_t      seqsource_header;
  multi1d<ForwardProp_t> forward_props;
};


//! Mega structure holding a full sequential prop (except gauge)
struct SequentialProp_t
{
  ChromaProp_t     seqprop_header;
  PropSink_t       sink_header;
  SeqSource_t      seqsource_header;
  multi1d<ForwardProp_t> forward_props;
};


//! Mega structure holding QQQ props (except gauge)
struct QQQBarcomp_t
{
  QQQBarcomp_t();                // default constructor
  bool             Dirac_basis;  // spin component basis
  multi1d<ForwardProp_t> forward_props;
};


//! Given a fermion action in string form, return the boundary
multi1d<int> getFermActBoundary(const std::string& fermact);


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

//! Multi Propagator header read
void read(XMLReader& xml, const std::string& path, ChromaMultiProp_t& header);


//! Propagator header writer
void write(XMLWriter& xml, const std::string& path, const ChromaProp_t& header);

//! Propagator header writer
void write(XMLWriter& xml, const std::string& path, const ChromaMultiProp_t& header);


//! SeqSource header read
void read(XMLReader& xml, const std::string& path, SeqSource_t& header);

//! SeqSource header writer
void write(XMLWriter& xml, const std::string& path, const SeqSource_t& header);


//! SeqProp header read
void read(XMLReader& xml, const std::string& path, ChromaSeqProp_t& header);

//! SeqProp header writer
void write(XMLWriter& xml, const std::string& path, const ChromaSeqProp_t& header);


//! ForwardProp reader
void read(XMLReader& xml, const std::string& path, ForwardProp_t& header);

//! ForwardProp writer
void write(XMLWriter& xml, const std::string& path, const ForwardProp_t& header);


//! SequentialSource reader
void read(XMLReader& xml, const std::string& path, SequentialSource_t& header);

//! SequentialSource writer
void write(XMLWriter& xml, const std::string& path, const SequentialSource_t& header);


//! SequentialProp reader
void read(XMLReader& xml, const std::string& path, SequentialProp_t& header);

//! SequentialProp writer
void write(XMLWriter& xml, const std::string& path, const SequentialProp_t& header);


//! QQQBarcomp reader
void read(XMLReader& xml, const std::string& path, QQQBarcomp_t& header);

//! QQQBarcomp writer
void write(XMLWriter& xml, const std::string& path, const QQQBarcomp_t& header);


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

// Write a Chroma Fermion Field (eg prop_component)
/*
 * \param file_xml     file header ( Read )
 * \param record_xml   xml holding propagator info ( Read )
 * \param fermion      fermion field( Read )
 * \param file         path ( Read )
 * \param volfmt       either QDPIO_SINGLEFILE, QDPIO_MULTIFILE ( Read )
 * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */    
void writeFermion(XMLBufferWriter& file_xml,
		  XMLBufferWriter& record_xml, const LatticeFermion& fermion,
		  const string& file, 
		  QDP_volfmt_t volfmt, QDP_serialparallel_t serpar);


// Read a Chroma Fermion Field
/*
 * \param file_xml     file header ( Write )
 * \param record_xml   xml holding propagator info ( Write )
 * \param fermion      The Fermion ( Write )
 * \param file         path ( Read )
 * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */    
void readFermion(XMLReader& file_xml,
		 XMLReader& record_xml, 
		 LatticeFermion& fermion,
		 const string& file, 
		 QDP_serialparallel_t serpar);

/*! @} */  // end of group io

}  // end namespace Chroma

#endif
