// -*- C++ -*-
// $Id: qprop_io.h,v 3.11 2007-11-27 23:01:26 kostas Exp $
/*! \file
 * \brief Routines associated with Chroma propagator IO
 */

#ifndef __qprop_io_h__
#define __qprop_io_h__

#include "io/xml_group_reader.h"
#include "io/enum_io/enum_qdpvolfmt_io.h"
#include "io/enum_io/enum_quarkspintype_io.h"

namespace Chroma 
{

  /*!
   * Chroma propagator support
   *
   * \ingroup io
   *
   * @{
   */


  //! Propagator source construction parameters
  struct PropSourceConst_t
  {
    PropSourceConst_t();              /*!< default constructor */

    multi1d<int> getTSrce() const;    /*!< return 4D coords of source (may not exist) */

    multi1d<int> getMom() const;    /*!< return the momentum of the source (may not exist) */

    GroupXML_t       source;          /*!< Holds source xml params*/

    int              j_decay;         /*!< decay direction */
    int              t_source;        /*!< source slice location */
  };

  //! Source-smearing parameters
  struct PropSourceSmear_t
  {
    PropSourceSmear_t();              /*!< default constructor */

    GroupXML_t       source;          /*!< Holds source xml params*/
    int              j_decay;         /*!< decay direction */
  };

  //! Sink-smearing parameters
  struct PropSinkSmear_t
  {
    PropSinkSmear_t();                /*!< default constructor */

    GroupXML_t       sink;            /*!< Holds sink xml params*/
    int              j_decay;         /*!< decay direction */
  };


  //! Multi-propagator parameters
  struct ChromaMultiProp_t 
  { 
    ChromaMultiProp_t();              /*!< default constructor */
    QuarkSpinType   quarkSpinType;    /*!< what spin components to compute */

    //! String holding XML of the FermionAction section
    GroupXML_t      fermact;          /*!< fermion action */
    GroupXML_t      invParam;         /*!< Inverter parameters */
 
    multi1d<Real>   MultiMasses;
  };


  //! Propagator parameters
  struct ChromaProp_t 
  { 
    ChromaProp_t();                   /*!< default constructor */
    QuarkSpinType   quarkSpinType;    /*!< why spin components to compute */

    // String holding XML of the FermionAction section
    GroupXML_t      fermact;          /*!< fermion action */
    bool            obsvP;            /*!< measure any observables (like Z_V, or mresP) on 5D prop */
  
    // String holding XML for auxiliary state information
    GroupXML_t      invParam;         /*!< Inverter parameters */
  };


  //! Sequential source parameters
  struct SeqSource_t
  {
    SeqSource_t();                    /*!< default constructor */

    GroupXML_t       seqsrc;          /*!< Sequential source xml */

    multi1d<int>     sink_mom;        /*!< sink momentum */
    int              t_sink;          /*!< time slice of sink */
    int              j_decay;         /*!< decay direction */
  };


  //! Mega structure holding a propagator source
  struct MakeSourceProp_t
  {
    PropSourceConst_t   source_header;
    std::string         gauge_header;
  };


  //! Mega structure holding a full forward prop (not sink smeared)
  struct Propagator_t
  {
    ChromaProp_t        prop_header;
    PropSourceConst_t   source_header;
    std::string         gauge_header;
  };


  //! Mega structure holding a full forward sink-smeared prop
  struct ForwardProp_t
  {
    PropSinkSmear_t     sink_header;
    ChromaProp_t        prop_header;
    PropSourceConst_t   source_header;
    std::string         gauge_header;
  };


  //! Mega structure holding a full sequential source
  struct SequentialSource_t
  {
    PropSinkSmear_t         sink_header;
    SeqSource_t             seqsource_header;
    multi1d<ForwardProp_t>  forward_props;
    std::string             gauge_header;
  };


  //! Mega structure holding a full sequential prop
  struct SequentialProp_t
  {
    ChromaProp_t           seqprop_header;
    PropSinkSmear_t        sink_header;
    SeqSource_t            seqsource_header;
    multi1d<ForwardProp_t> forward_props;
    std::string            gauge_header;
  };


  //! Mega structure holding a full sequential prop that is source smeared
  struct SeqPropSourceSmeared_t
  {
    PropSourceSmear_t      smeared_seqprop_header;
    ChromaProp_t           seqprop_header;
    PropSinkSmear_t        sink_header;
    SeqSource_t            seqsource_header;
    multi1d<ForwardProp_t> forward_props;
    std::string            gauge_header;
  };


  //! Hold source and sink spin indices for a sparse QQQ file
  struct QQQSpinIndices_t
  {
    multi1d<int>  source;
    multi1d<int>  sink;
  };


  //! Mega structure holding QQ diquark object
  struct QQDiquark_t
  {
    QQDiquark_t();                 /*!< default constructor */
    bool             Dirac_basis;  /*!< spin component basis */
    multi1d<ForwardProp_t> forward_props;
  };


  //! Mega structure holding QQQ props
  struct QQQBarcomp_t
  {
    QQQBarcomp_t();                // default constructor
    multi1d<QQQSpinIndices_t> spin_indices;  // spin indices
    bool             sparseP;      // Is this data sparsely stored?
    bool             Dirac_basis;  // spin component basis
    multi1d<ForwardProp_t> forward_props;
  };


  //! Mega structure holding QQbar props
  struct QQbarMescomp_t
  {
    QQbarMescomp_t();                // default constructor
    bool             Dirac_basis;    // spin component basis
    multi1d<ForwardProp_t> forward_props;
  };


  //! Given a fermion action in string form, return the Mass
  Real getMass(const GroupXML_t& fermact);

  //! Given a fermion action in string form, return the boundary
  multi1d<int> getFermActBoundary(const GroupXML_t& fermact);

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


  //! MakeSourceProp reader
  void read(XMLReader& xml, const std::string& path, MakeSourceProp_t& header);

  //! MakeSourceProp writer
  void write(XMLWriter& xml, const std::string& path, const MakeSourceProp_t& header);


  //! Propagator_t reader
  void read(XMLReader& xml, const std::string& path, Propagator_t& header);

  //! Propagator_t writer
  void write(XMLWriter& xml, const std::string& path, const Propagator_t& header);


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


  //! Source/sink spin indices
  void read(XMLReader& xml, const string& path, QQQSpinIndices_t& input);

  //! Source/sink spin indices
  void write(XMLWriter& xml, const string& path, const QQQSpinIndices_t& input);


  //! QQDiquark reader
  void read(XMLReader& xml, const std::string& path, QQDiquark_t& header);

  //! QQDiquark writer
  void write(XMLWriter& xml, const std::string& path, const QQDiquark_t& header);


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
