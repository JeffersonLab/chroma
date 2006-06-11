// -*- C++ -*-
// $Id: qprop_io.h,v 3.1 2006-06-11 06:30:33 edwards Exp $
/*! \file
 * \brief Routines associated with Chroma propagator IO
 */

#ifndef __qprop_io_h__
#define __qprop_io_h__

#include "io/enum_io/enum_invtype_io.h"
#include "io/enum_io/enum_qdpvolfmt_io.h"
#include "io/enum_io/enum_quarkspintype_io.h"

namespace Chroma {

/*!
 * Chroma propagator support
 *
 * \ingroup io
 *
 * @{
 */


//! Propagator source construction header
struct PropSourceConst_t
{
  PropSourceConst_t();              /*!< default constructor */

  multi1d<int> getTSrce();          /*!< return 4D coords of source (may not exist) */

  std::string      source;          /*!< string holding source xml */
  std::string      source_type;     /*!< source type */

  int              j_decay;         /*!< decay direction */
  int              t_source;        /*!< source slice location */
};

//! Propagator source smearing header
struct PropSourceSmear_t
{
  PropSourceSmear_t();              /*!< default constructor */

  std::string      source;          /*!< string holding source xml */
  std::string      source_type;     /*!< source type */

  int              j_decay;         /*!< decay direction */
};

//! Propagator sink header
struct PropSinkSmear_t
{
  PropSinkSmear_t();                /*!< default constructor */

  std::string      sink;            /*!< string holding sink smearing xml */
  std::string      sink_type;       /*!< sink type */

  int              j_decay;         /*!< decay direction */
};


//! Multiple propagator header
struct ChromaMultiProp_t 
{ 
  ChromaMultiProp_t();              /*!< default constructor */
  QuarkSpinType   quarkSpinType;    /*!< what spin components to compute */

  //! String holding XML of the FermionAction section
  std::string     fermact;          /*!< fermion action */

  MultiInvertParam_t   invParam;    /*!< Inverter parameters */
 
  multi1d<Real>   MultiMasses;
};


//! Propagator header with action
struct ChromaProp_t 
{ 
  ChromaProp_t();                   /*!< default constructor */
  QuarkSpinType   quarkSpinType;    /*!< why spin components to compute */

  // String holding XML of the FermionAction section
  std::string     fermact;          /*!< fermion action */

  bool            obsvP;            /*!< measure any observables (like Z_V, or mresP) on 5D prop */
  int             numRetries;       /*!< number of calls to qprop for each source component */
  
  // String holding XML for auxiliary state information
  InvertParam_t   invParam;         /*!< Inverter parameters */
};


//! Structure for writing to seqsource files
struct SeqSource_t
{
  SeqSource_t();                    /*!< default constructor */

  std::string      seqsrc;          /*!< string holding sequential source xml */
  std::string      seqsrc_type;     /*!< sequential source type */

  multi1d<int>     sink_mom;        /*!< sink momentum */
  int              t_sink;          /*!< time slice of sink */
  int              j_decay;         /*!< decay direction */
};


//! Structure for writing to seqprop files
struct ChromaSeqProp_t
{
  QuarkSpinType    quarkSpinType;   // which spin components to compute
  InvertParam_t    invParam;
  std::string      seq_src;
  multi1d<int>     sink_mom;
  int              t_sink;
};


//! Mega structure holding a full forward prop (except gauge)
struct ForwardProp_t
{
  PropSinkSmear_t     sink_header;
  ChromaProp_t        prop_header;
  PropSourceConst_t   source_header;
};


//! Mega structure holding a full sequential source (except gauge)
struct SequentialSource_t
{
  PropSinkSmear_t         sink_header;
  SeqSource_t             seqsource_header;
  multi1d<ForwardProp_t>  forward_props;
};


//! Mega structure holding a full sequential prop (except gauge)
struct SequentialProp_t
{
  ChromaProp_t           seqprop_header;
  PropSinkSmear_t        sink_header;
  SeqSource_t            seqsource_header;
  multi1d<ForwardProp_t> forward_props;
};


//! Mega structure holding QQQ props (except gauge)
struct QQQBarcomp_t
{
  QQQBarcomp_t();                // default constructor
  bool             Dirac_basis;  // spin component basis
  multi1d<ForwardProp_t> forward_props;
};


//! Mega structure holding QQbar props (except gauge)
struct QQbarMescomp_t
{
  QQbarMescomp_t();                // default constructor
  bool             Dirac_basis;    // spin component basis
  multi1d<ForwardProp_t> forward_props;
};


//! Given a fermion action in string form, return the boundary
multi1d<int> getFermActBoundary(const std::string& fermact);

//! Propagator source read
void read(XMLReader& xml, const std::string& path, PropSourceConst_t& header);

//! Propagator source writer
void write(XMLWriter& xml, const std::string& path, const PropSourceConst_t& header);


//! Propagator source smearing read
void read(XMLReader& xml, const std::string& path, PropSourceSmear_t& header);

//! Propagator source smearing writer
void write(XMLWriter& xml, const std::string& path, const PropSourceSmear_t& header);


//! Propagator sink reader
void read(XMLReader& xml, const std::string& path, PropSinkSmear_t& header);

//! Propagator sink writer
void write(XMLWriter& xml, const std::string& path, const PropSinkSmear_t& header);


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


//! QQbarMescomp reader
void read(XMLReader& xml, const std::string& path, QQbarMescomp_t& header);

//! QQbarMescomp writer
void write(XMLWriter& xml, const std::string& path, const QQbarMescomp_t& header);


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
