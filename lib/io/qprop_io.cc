// $Id: qprop_io.cc,v 1.12 2004-04-01 18:09:58 edwards Exp $
/*! \file
 * \brief Routines associated with Chroma propagator IO
 */

#include "chromabase.h"
#include "io/param_io.h"
#include "io/qprop_io.h"

#include <string>
using std::string;


// Initialize header with default values
void initHeader(PropSource_t& header)
{
  header.version     = 5;
  header.source_type = SRC_TYPE_POINT_SOURCE;
  header.wave_state  = WAVE_TYPE_S_WAVE;

  initHeader(header.sourceSmearParam);

  header.j_decay         = 0;
  header.direction       = 0;
  header.laplace_power   = 0;
  header.disp_length     = 0;
  header.disp_dir        = 0;
  header.link_smear_fact = 0;
  header.link_smear_num  = 0;
  header.t_source.resize(Nd);
  header.t_source        = 0;
  header.nrow            = Layout::lattSize();
}

// Initialize header with default values
void initHeader(PropSink_t& header)
{
  header.version       = 4;
  header.sink_type     = SNK_TYPE_POINT_SINK;
  header.wave_state    = WAVE_TYPE_S_WAVE;

  initHeader(header.sinkSmearParam);

  header.direction       = 0;
  header.laplace_power   = 0;
  header.disp_length     = 0;
  header.disp_dir        = 0;
  header.link_smear_fact = 0;
  header.link_smear_num  = 0;
  header.nrow            = Layout::lattSize();
}

// Initialize header with default values
void initHeader(PropSink_t& header, const PropSource_t& source)
{
  header.version = 4;

  // Convert the source to a sink type
  switch(source.source_type)
  {
  case SRC_TYPE_POINT_SOURCE:
    header.sink_type = SNK_TYPE_POINT_SINK;
    break;
  case SRC_TYPE_WALL_SOURCE:
    header.sink_type = SNK_TYPE_WALL_SINK;
    break;
  case SRC_TYPE_SHELL_SOURCE:
    header.sink_type = SNK_TYPE_SHELL_SINK;
    break;
  case SRC_TYPE_BNDST_SOURCE:
    header.sink_type = SNK_TYPE_BNDST_SINK;
    break;
  case SRC_TYPE_POINT_AND_BNDST_SOURCE:
    header.sink_type = SNK_TYPE_POINT_AND_BNDST_SINK;
    break;
  case SRC_TYPE_SHELL_AND_BNDST_SOURCE:
    header.sink_type = SNK_TYPE_SHELL_AND_BNDST_SINK;
    break;
  case SRC_TYPE_POINT_AND_SHELL_AND_BNDST_SOURCE:
    header.sink_type = SNK_TYPE_POINT_AND_SHELL_AND_BNDST_SINK;
    break;
  }

  // description (same order):
  // wave_state: "S_WAVE", "P_WAVE", or "D_WAVE"
  // smearing params
  //   smearing width
  //   number of iteration for smearing
  //   number of iteration for smearing
  // decay direction
  // power of laplacian operator
  // displacement length
  // displacement direction: x(0),y(1),z(2)
  // sink_dir: direction of derivative at sink
  // smearing factor
  // number of smearing hits
  // Maximum number of blocking/smearing iterations
  // Blocking/smearing accuracy

  header.wave_state      = source.wave_state;
  header.sinkSmearParam  = source.sourceSmearParam;
  header.direction       = source.direction;
  header.laplace_power   = source.laplace_power;
  header.disp_length     = source.disp_length;
  header.disp_dir        = source.disp_dir;
  header.link_smear_fact = source.link_smear_fact;
  header.link_smear_num  = source.link_smear_num;
  header.nrow            = Layout::lattSize();
}


// Initialize header with default values
void initHeader(ChromaProp_t& header)
{
  header.version     = 4;
  header.nrow        = Layout::lattSize();

  initHeader(header.anisoParam);
  initHeader(header.chiralParam);
}



// Source header read
void read(XMLReader& xml, const string& path, PropSource_t& header)
{
  XMLReader paramtop(xml, path);

  initHeader(header);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
    /**************************************************************************/
  case 5:
    /**************************************************************************/
    break;

  default:
    /**************************************************************************/
    QDPIO::cerr << "PropSource parameter version " << version 
		<< " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "source_type", header.source_type);
  read(paramtop, "wave_state", header.wave_state);
  read(paramtop, "j_decay",  header.j_decay);
  read(paramtop, "direction",  header.direction);
  read(paramtop, "t_source", header.t_source);

  if (header.source_type == SRC_TYPE_SHELL_SOURCE)
  {
    XMLReader shelltop(paramtop, "ShellSource");

    read(shelltop, "SourceSmearingParam", header.sourceSmearParam);
    read(shelltop, "laplace_power", header.laplace_power);
    read(shelltop, "link_smear_fact", header.link_smear_fact);
    read(shelltop, "link_smear_num", header.link_smear_num);
    read(shelltop, "disp_length", header.disp_length);
    read(shelltop, "disp_dir", header.disp_dir);
  }		
  read(paramtop, "nrow", header.nrow);
}

// Source header read
void read(XMLReader& xml, const string& path, PropSink_t& header)
{
  XMLReader paramtop(xml, path);

  initHeader(header);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
    /**************************************************************************/
  case 4:
    /**************************************************************************/
    break;

  default:
    /**************************************************************************/
    QDPIO::cerr << "PropSink parameter version " << version 
		<< " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "sink_type", header.sink_type);
  read(paramtop, "wave_state", header.wave_state);
  read(paramtop, "direction",  header.direction);

  if (header.sink_type == SNK_TYPE_SHELL_SINK)
  {
    XMLReader shelltop(paramtop, "ShellSink");

    read(shelltop, "SinkSmearingParam", header.sinkSmearParam);
    read(shelltop, "laplace_power", header.laplace_power);
    read(shelltop, "link_smear_fact", header.link_smear_fact);
    read(shelltop, "link_smear_num", header.link_smear_num);
    read(shelltop, "disp_length", header.disp_length);
    read(shelltop, "disp_dir", header.disp_dir);
  }		
  read(paramtop, "nrow", header.nrow);
}

//! SeqPropagator header reader
void read(XMLReader& xml, const string& path, ChromaSeqProp_t& param)
{
  XMLReader paramtop(xml, path);

//  initHeader(param);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
    /**************************************************************************/
  case 1:
    param.nonRelSeqProp = false;
    /**************************************************************************/
    break;

    /**************************************************************************/
  case 2:
    read(paramtop, "nonRelSeqProp", param.nonRelSeqProp);
    break;

  default:
    /**************************************************************************/
    QDPIO::cerr << "ChromaSeqProp parameter version " << version 
		<< " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "Seq_src", param.Seq_src);
  read(paramtop, "InvertParam", param.invParam);
  read(paramtop, "t_sink", param.t_sink);
  read(paramtop, "sink_mom", param.sink_mom);
  read(paramtop, "nrow", param.nrow);
}


// Forward propagator header read
void read(XMLReader& xml, const string& path, ChromaProp_t& param)
{
  XMLReader paramtop(xml, path);

  initHeader(param);

  int version;
  read(paramtop, "version", version);

  switch (version) 
  {
    /**************************************************************************/
  case 4:
    /**************************************************************************/
    break;

  default:
    /**************************************************************************/
    QDPIO::cerr << "ChromaProp parameter version " << version 
		<< " unsupported." << endl;
    QDP_abort(1);
  }

  read(paramtop, "FermTypeP", param.FermTypeP);
  read(paramtop, "FermAct", param.FermAct);

  if (paramtop.count("Mass") != 0)
  {
    read(paramtop, "Mass", param.Mass);

    if (paramtop.count("Kappa") != 0)
    {
      QDPIO::cerr << "Error: found both a Kappa and a Mass tag" << endl;
      QDP_abort(1);
    }
  }
  else if (paramtop.count("Kappa") != 0)
  {
    Real Kappa;
    read(paramtop, "Kappa", Kappa);

    param.Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
  }
  else
  {
    QDPIO::cerr << "Error: neither Mass or Kappa found" << endl;
    QDP_abort(1);
  }    

  if (paramtop.count("AnisoParam") != 0)
    read(paramtop, "AnisoParam", param.anisoParam);

  if (paramtop.count("ChiralParam") != 0)
    read(paramtop, "ChiralParam", param.chiralParam);

  read(paramtop, "InvertParam", param.invParam);
  read(paramtop, "boundary", param.boundary);
  read(paramtop, "nrow", param.nrow);
}



// Source header writer
void write(XMLWriter& xml, const string& path, const PropSource_t& header)
{
  push(xml, path);

  write(xml, "version", header.version);
  write(xml, "source_type", header.source_type);
  write(xml, "wave_state", header.wave_state);
  write(xml, "j_decay",  header.j_decay);
  write(xml, "direction",  header.direction);
  write(xml, "t_source",  header.t_source);

  if (header.source_type == SRC_TYPE_SHELL_SOURCE)
  {
    push(xml, "ShellSource");
    write(xml, "SourceSmearingParam", header.sourceSmearParam);
    write(xml, "laplace_power", header.laplace_power);
    write(xml, "link_smear_fact", header.link_smear_fact);
    write(xml, "link_smear_num", header.link_smear_num);
    write(xml, "disp_length", header.disp_length);
    write(xml, "disp_dir", header.disp_dir);
    pop(xml);
  }

  write(xml, "nrow",  header.nrow);

  pop(xml);
}


// Source header writer
void write(XMLWriter& xml, const string& path, const PropSink_t& header)
{
  push(xml, path);

  write(xml, "version", header.version);
  write(xml, "sink_type", header.sink_type);
  write(xml, "wave_state", header.wave_state);
  write(xml, "direction",  header.direction);

  if (header.sink_type == SNK_TYPE_SHELL_SINK)
  {
    push(xml, "ShellSink");
    write(xml, "SinkSmearingParam", header.sinkSmearParam);
    write(xml, "laplace_power", header.laplace_power);
    write(xml, "link_smear_fact", header.link_smear_fact);
    write(xml, "link_smear_num", header.link_smear_num);
    write(xml, "disp_length", header.disp_length);
    write(xml, "disp_dir", header.disp_dir);
    pop(xml);
  }

  write(xml, "nrow",  header.nrow);

  pop(xml);
}


// Write propagator inversion parameters
void write(XMLWriter& xml, const string& path, const ChromaProp_t& header)
{
  push(xml, path);

  write(xml, "version", header.version);
  write(xml, "FermTypeP", header.FermTypeP);
  write(xml, "FermAct", header.FermAct);
  write(xml, "Mass", header.Mass);
  write(xml, "AnisoParam", header.anisoParam);
  write(xml, "ChiralParam", header.chiralParam);
  write(xml, "InvertParam", header.invParam);
  write(xml, "boundary", header.boundary);
  write(xml, "nrow", header.nrow);

  pop(xml);
}


//! SeqPropagator header writer
void write(XMLWriter& xml, const string& path, const ChromaSeqProp_t& param)
{
  push(xml, path);

  int version = 2;
  write(xml, "version", version);
  write(xml, "nonRelSeqProp", param.nonRelSeqProp);
  write(xml, "Seq_src", param.Seq_src);
  write(xml, "InvertParam", param.invParam);
  write(xml, "t_sink", param.t_sink);
  write(xml, "sink_mom", param.sink_mom);
  write(xml, "nrow", param.nrow);

  pop(xml);
}



// Write a Chroma propagator
/*
 * \param file_xml     file header ( Read )
 * \param record_xml   xml holding propagator info ( Read )
 * \param quark_prop   propagator ( Read )
 * \param file         path ( Read )
 * \param volfmt       either QDPIO_SINGLEFILE, QDPIO_MULTIFILE ( Read )
 * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */    
void writeQprop(XMLBufferWriter& file_xml,
		XMLBufferWriter& record_xml, const LatticePropagator& quark_prop,
		const string& file, 
		QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
{
  QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
  write(to,record_xml,quark_prop);
  close(to);
}

// Write a Chroma propagator
/*
 * \param file_xml     file header ( Read )
 * \param header       structure holding propagator info ( Read )
 * \param quark_prop   propagator ( Read )
 * \param file         path ( Read )
 * \param volfmt       either QDPIO_SINGLEFILE, QDPIO_MULTIFILE ( Read )
 * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */    
void writeQprop(XMLBufferWriter& file_xml,
		const ChromaProp_t& header, const LatticePropagator& quark_prop,
		const string& file, 
		QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
{
  XMLBufferWriter record_xml;
  write(record_xml, "Propagator", header);
  writeQprop(file_xml, record_xml, quark_prop, file, volfmt, serpar);
}




// Read a Chroma propagator
/*
 * \param file_xml     file header ( Write )
 * \param record_xml   xml holding propagator info ( Write )
 * \param quark_prop   propagator ( Write )
 * \param file         path ( Read )
 * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */    
void readQprop(XMLReader& file_xml,
	       XMLReader& record_xml, LatticePropagator& quark_prop,
	       const string& file, 
	       QDP_serialparallel_t serpar)
{
  QDPFileReader to(file_xml,file,serpar);
  read(to,record_xml,quark_prop);
  close(to);
}

// Read a Chroma propagator
/*
 * \param file_xml     file header ( Write )
 * \param header       structure holding propagator info ( Write )
 * \param quark_prop   propagator ( Write )
 * \param file         path ( Read )
 * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
 */    
void readQprop(XMLReader& file_xml,
	       ChromaProp_t& header, LatticePropagator& quark_prop,
	       const string& file, 
	       QDP_serialparallel_t serpar)
{
  XMLReader record_xml;
  readQprop(file_xml, record_xml, quark_prop, file, serpar);
  read(record_xml, "/Propagator", header);  // extract header from xml
}




